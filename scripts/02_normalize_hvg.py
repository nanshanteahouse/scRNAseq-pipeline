#!/usr/bin/env python3
"""
Step 02: 归一化 + HVG 选择
=============================
关键顺序 (继承 GSE169109):
  1. raw counts → 找 HVG (seurat_v3, batch-aware)
  2. normalize_total → log1p
  3. adata.raw = adata  (保留全基因)
  4. X 只保留 HVG 用于下游降维

输入: 02_qc.h5ad
输出: 02_normalized.h5ad (X = log1p(normalized) on HVGs, .raw = 全基因)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write
import scanpy as sc

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("02_norm_hvg", os.path.join(CFG.log_dir, "02_normalize_hvg.log"))
    log.info("Step 02: 归一化 + HVG 选择")

    adata = sc.read(CFG.qc_h5ad)
    log.info("加载: %d 细胞 × %d 基因", adata.n_obs, adata.n_vars)

    # Step A: 在 raw counts 上找 HVG (批次感知)
    log.info("选择 top %d HVGs (flavor=%s, batch_key=%s)...",
             CFG.n_top_genes, CFG.hvg_flavor, CFG.hvg_batch_key)
    batch_key = CFG.hvg_batch_key if CFG.hvg_batch_key in adata.obs else None
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=CFG.n_top_genes,
        flavor=CFG.hvg_flavor,
        batch_key=batch_key,
        inplace=True,
    )
    n_hvg = adata.var['highly_variable'].sum()
    log.info("HVG 数量: %d", n_hvg)

    # Step B: 归一化
    log.info("归一化 (target_sum=%.0f) + log1p...", CFG.normalize_target_sum)
    sc.pp.normalize_total(adata, target_sum=CFG.normalize_target_sum)
    sc.pp.log1p(adata)

    # Step C: 保留全基因到 .raw
    adata.raw = adata
    log.info(".raw 已保存 (全基因: %d vars)", adata.raw.n_vars)

    # Step D: X 只保留 HVG
    adata = adata[:, adata.var['highly_variable']].copy()
    log.info("X 缩小到 HVGs: %s", adata.shape)

    safe_write(adata, CFG.norm_h5ad)
    log.info("Step 02 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
