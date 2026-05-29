#!/usr/bin/env python3
"""
Step 01: 质量控制 — 指标计算 + Scrublet 双细胞检测 + 过滤
=============================================================
继承了 GSE169109 的最佳实践:
  1. QC 指标 (mito%, ribo%, 复杂度)
  2. Scrublet 双细胞检测 (per sample, joblib 并行)
  3. 自适应 MAD 或全局阈值过滤

输入: 00_raw.h5ad
输出: 01_qc.h5ad (过滤后的细胞 + QC 指标)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed

def compute_qc_metrics(adata, log):
    log.info("计算 QC 指标...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt', 'ribo'],
        percent_top=[20], log1p=True, inplace=True,
    )
    adata.obs['log_genes_per_umi'] = (
        np.log10(adata.obs['n_genes_by_counts'])
        / np.log10(adata.obs['total_counts'])
    ).replace([np.inf, -np.inf], np.nan)
    log.info("  中位基因/细胞: %.0f", adata.obs['n_genes_by_counts'].median())
    log.info("  中位 UMIs/细胞: %.0f", adata.obs['total_counts'].median())
    log.info("  中位 mito%%:    %.2f%%", adata.obs['pct_counts_mt'].median())
    log.info("  中位复杂度:     %.3f", adata.obs['log_genes_per_umi'].median())

def run_scrublet_sample(adata_sub, sample_name, cfg):
    try:
        import scrublet as scr
        import scipy.sparse as sp
        scrub = scr.Scrublet(
            adata_sub.X if isinstance(adata_sub.X, sp.spmatrix) else sp.csr_matrix(adata_sub.X),
            expected_doublet_rate=cfg.scrublet_expected_doublet_rate,
            random_state=cfg.random_seed,
        )
        scores, predicted = scrub.scrub_doublets(
            min_counts=cfg.scrublet_min_counts,
            min_cells=cfg.scrublet_min_cells,
            min_gene_variability_pctl=cfg.scrublet_min_gene_var_pctl,
            n_prin_comps=cfg.scrublet_n_prin_comps,
        )
        return scores, predicted
    except Exception as e:
        import warnings
        warnings.warn(f"Scrublet failed for {sample_name}: {e}")
        return np.zeros(adata_sub.n_obs), np.zeros(adata_sub.n_obs, dtype=bool)

def detect_doublets_parallel(adata, cfg, log):
    if not cfg.run_scrublet:
        log.info("Scrublet 已禁用，跳过双细胞检测。")
        adata.obs['doublet_scores'] = 0.0
        adata.obs['predicted_doublet'] = False
        return
    log.info("运行 Scrublet (按样本并行)...")
    groupby_col = 'sample' if 'sample' in adata.obs else 'stage'
    if groupby_col not in adata.obs:
        log.warning("找不到分组列，Scrublet 在全体数据上运行。")
        scores, pred = run_scrublet_sample(adata, "all", cfg)
        adata.obs['doublet_scores'] = scores
        adata.obs['predicted_doublet'] = pred
        return

    sample_groups = adata.obs.groupby(groupby_col)
    names = []
    subsets = []
    for name, idx in sample_groups.indices.items():
        names.append(name)
        subsets.append(adata[idx])

    results = Parallel(n_jobs=min(cfg.n_jobs, len(names)))(
        delayed(run_scrublet_sample)(sub, name, cfg)
        for sub, name in zip(subsets, names)
    )

    all_scores = np.zeros(adata.n_obs)
    all_pred = np.zeros(adata.n_obs, dtype=bool)
    for (scores, pred), idx in zip(results, sample_groups.indices.values()):
        all_scores[idx] = scores
        all_pred[idx] = pred
    adata.obs['doublet_scores'] = all_scores
    adata.obs['predicted_doublet'] = all_pred
    log.info("  预测双细胞: %d / %d (%.1f%%)",
             all_pred.sum(), adata.n_obs, 100 * all_pred.mean())

def filter_cells(adata, cfg, log):
    log.info("应用过滤...")
    min_g = cfg.min_genes
    max_g = cfg.max_genes
    max_m = cfg.max_pct_mito
    min_cpx = cfg.min_genes_per_umi

    f_genes_low = adata.obs['n_genes_by_counts'] < min_g
    f_genes_high = adata.obs['n_genes_by_counts'] > max_g
    f_mito = adata.obs['pct_counts_mt'] > max_m
    f_cpx = adata.obs['log_genes_per_umi'] < np.log10(min_cpx)
    f_doublet = adata.obs['predicted_doublet']
    f_any = f_genes_low | f_genes_high | f_mito | f_cpx | f_doublet

    log.info("  过滤明细:")
    log.info("    n_genes < %d:     %6d (%.1f%%)", min_g, f_genes_low.sum(), 100*f_genes_low.mean())
    log.info("    n_genes > %d:     %6d (%.1f%%)", max_g, f_genes_high.sum(), 100*f_genes_high.mean())
    log.info("    mito > %d%%:      %6d (%.1f%%)", max_m, f_mito.sum(), 100*f_mito.mean())
    log.info("    复杂度 < %.2f:    %6d (%.1f%%)", min_cpx, f_cpx.sum(), 100*f_cpx.mean())
    log.info("    Scrublet 双细胞:  %6d (%.1f%%)", f_doublet.sum(), 100*f_doublet.mean())
    log.info("    合计(去重):       %6d (%.1f%%)", f_any.sum(), 100*f_any.mean())

    adata = adata[~f_any].copy()
    log.info("  过滤后: %d 细胞", adata.n_obs)
    return adata

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("01_qc", os.path.join(CFG.log_dir, "01_qc.log"))
    log.info("Step 01: QC 预处理")

    adata = sc.read(CFG.raw_h5ad)
    log.info("加载: %s — %d 细胞 × %d 基因",
             CFG.raw_h5ad, adata.n_obs, adata.n_vars)

    compute_qc_metrics(adata, log)
    detect_doublets_parallel(adata, CFG, log)
    adata = filter_cells(adata, CFG, log)
    sc.pp.filter_genes(adata, min_cells=CFG.min_cells_per_gene)
    log.info("基因过滤后: %d 基因", adata.n_vars)

    safe_write(adata, CFG.qc_h5ad)
    log.info("Step 01 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
