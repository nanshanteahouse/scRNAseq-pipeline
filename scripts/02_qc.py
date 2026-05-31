#!/usr/bin/env python3
"""
Step 02: QC filtering (doublets already removed in Step 01)
============================================================
继承了 GSE169109 的最佳实践:
  1. QC 指标 (mito%, ribo%, 复杂度)
  2. 过滤 predicted_doublet 细胞 (由 Step 01 产生)
  3. 自适应 MAD 或全局阈值过滤

输入: 01_doublet.h5ad (含 doublet_scores, predicted_doublet 列)
输出: 02_qc.h5ad (过滤后的细胞 + QC 指标)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write
import numpy as np
import pandas as pd
import scanpy as sc

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



def filter_cells(adata, cfg, log):
    # Step 1: Filter predicted doublets first (from Step 01)
    n_before = adata.n_obs
    adata = adata[~adata.obs['predicted_doublet']].copy()
    n_doublet = n_before - adata.n_obs
    log.info("双细胞过滤: 移除 %d 个预测双细胞 (%.1f%%)",
             n_doublet, 100 * n_doublet / n_before if n_before else 0)

    # Step 2: QC threshold filtering
    log.info("应用 QC 过滤...")
    min_g = cfg.min_genes
    max_g = cfg.max_genes
    max_m = cfg.max_pct_mito
    min_cpx = cfg.min_genes_per_umi

    f_genes_low  = adata.obs['n_genes_by_counts'] < min_g
    f_genes_high = adata.obs['n_genes_by_counts'] > max_g
    f_mito       = adata.obs['pct_counts_mt'] > max_m
    f_cpx        = adata.obs['log_genes_per_umi'] < np.log10(min_cpx)
    f_any = f_genes_low | f_genes_high | f_mito | f_cpx

    log.info("  过滤明细:")
    log.info("    n_genes < %d:     %6d (%.1f%%)", min_g, f_genes_low.sum(), 100*f_genes_low.mean())
    log.info("    n_genes > %d:     %6d (%.1f%%)", max_g, f_genes_high.sum(), 100*f_genes_high.mean())
    log.info("    mito > %d%%:      %6d (%.1f%%)", max_m, f_mito.sum(), 100*f_mito.mean())
    log.info("    复杂度 < %.2f:    %6d (%.1f%%)", min_cpx, f_cpx.sum(), 100*f_cpx.mean())
    log.info("    合计(去重):       %6d (%.1f%%)", f_any.sum(), 100*f_any.mean())

    adata = adata[~f_any].copy()
    log.info("  QC 过滤后: %d 细胞", adata.n_obs)
    return adata

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("02_qc", os.path.join(CFG.log_dir, "02_qc.log"))
    log.info("Step 02: QC 过滤 (双细胞已在 Step 01 中移除)")

    input_path = os.path.join(CFG.h5ad_dir, "01_doublet.h5ad")
    adata = sc.read(input_path)
    log.info("加载: %s — %d 细胞 × %d 基因",
             input_path, adata.n_obs, adata.n_vars)

    compute_qc_metrics(adata, log)
    adata = filter_cells(adata, CFG, log)
    sc.pp.filter_genes(adata, min_cells=CFG.min_cells_per_gene)
    log.info("基因过滤后: %d 基因", adata.n_vars)

    safe_write(adata, CFG.qc_h5ad)
    log.info("Step 02 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
