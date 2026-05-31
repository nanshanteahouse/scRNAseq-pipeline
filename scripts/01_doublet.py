#!/usr/bin/env python3
"""
Step 01a: Scrublet 双细胞检测 (per sample, joblib 并行)
=========================================================
从 01_qc.py 中独立出的 Scrublet 步骤，不含 QC 指标计算或过滤。

输入: 00_raw.h5ad
输出: 01_doublet.h5ad (含 doublet_scores / predicted_doublet 列)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write
import numpy as np
import scanpy as sc
from joblib import Parallel, delayed


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
        log.warning("找不到分组列(%s)，Scrublet 在全体数据上运行。", groupby_col)
        scores, pred = run_scrublet_sample(adata, "all", cfg)
        adata.obs['doublet_scores'] = scores
        adata.obs['predicted_doublet'] = pred
        log.info("  预测双细胞: %d / %d (%.1f%%)",
                 pred.sum(), adata.n_obs, 100 * pred.mean())
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
    for (scores, pred), name, idx in zip(results, names, sample_groups.indices.values()):
        all_scores[idx] = scores
        all_pred[idx] = pred
        log.info("  样本 %s: %d / %d 双细胞 (%.1f%%)",
                 name, pred.sum(), len(idx), 100 * pred.mean())

    adata.obs['doublet_scores'] = all_scores
    adata.obs['predicted_doublet'] = all_pred
    log.info("  总计预测双细胞: %d / %d (%.1f%%)",
             all_pred.sum(), adata.n_obs, 100 * all_pred.mean())


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("01_doublet", os.path.join(CFG.log_dir, "01_doublet.log"))
    log.info("Step 01a: Scrublet 双细胞检测")

    adata = sc.read(CFG.raw_h5ad)
    log.info("加载: %s — %d 细胞 × %d 基因",
             CFG.raw_h5ad, adata.n_obs, adata.n_vars)

    detect_doublets_parallel(adata, CFG, log)

    out_path = os.path.join(CFG.h5ad_dir, "01_doublet.h5ad")
    safe_write(adata, out_path)
    log.info("Step 01a 完成, 耗时 %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
