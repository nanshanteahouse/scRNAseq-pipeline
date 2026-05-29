#!/usr/bin/env python3
"""
Step 00: 加载原始 scRNA-seq 数据
===================================
支持三种输入格式:
  1. 10X MTX (CellRanger 输出): sc.read_10x_mtx()
  2. CSV 矩阵 + 元数据文件:     mmread() + pandas
  3. 已有 h5ad:                sc.read()

输出: 00_raw.h5ad
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("00_load", os.path.join(CFG.log_dir, "00_load.log"))
    log.info("Step 00: 加载原始数据")
    log.info("格式: %s", CFG.data_format)

    if os.path.exists(CFG.raw_h5ad):
        log.info("跳过: %s 已存在。删除它以强制重新加载。", CFG.raw_h5ad)
        return

    # ── 3 种加载方式 ──────────────────────────────────────────────
    if CFG.data_format == "10X_mtx":
        log.info("从 MTX 加载 (前缀='%s') ...", CFG.mtx_prefix)
        adata = sc.read_10x_mtx(
            CFG.mtx_dir,
            var_names='gene_symbols',
            prefix=CFG.mtx_prefix,
            cache=True,
            gex_only=False,
        )
        log.info("加载完成: %d 细胞 × %d 基因", adata.n_obs, adata.n_vars)

        # 解析 barcode 后缀 → 样本/阶段映射
        if CFG.has_sample_mapping() or CFG.has_stage_mapping():
            bc_suffix = (
                adata.obs_names.to_series()
                .str.extract(r'-(\d+)$')[0]
                .astype(int)
            )
            if CFG.has_sample_mapping():
                adata.obs['sample'] = bc_suffix.map(CFG.sample_map).values
            if CFG.has_stage_mapping():
                adata.obs['stage'] = bc_suffix.map(CFG.stage_map).values
                if CFG.stage_order:
                    adata.obs['stage'] = pd.Categorical(
                        adata.obs['stage'],
                        categories=CFG.stage_order,
                        ordered=True,
                    )
            log.info("样本映射已应用。样本分布:")
            if 'sample' in adata.obs:
                for s, cnt in adata.obs['sample'].value_counts().items():
                    log.info("  %-20s %5d cells", s, cnt)

        # 清理 gene_ids 列（如果有）
        if 'gene_ids' in adata.var:
            adata.var.drop(columns=['gene_ids'], inplace=True)

    elif CFG.data_format == "csv_matrix":
        log.info("从 CSV 矩阵加载: %s", CFG.matrix_file)
        mtx = mmread(CFG.matrix_file)
        log.info("矩阵形状: %s, nnz=%d", mtx.shape, mtx.nnz)
        mtx.data = mtx.data.astype(np.float32)
        mtx = mtx.T.tocsr()

        genes = pd.read_csv(CFG.features_file)
        gene_names = genes.iloc[:, 0].values.astype(str)
        gene_names = pd.Index(gene_names)
        if gene_names.duplicated().any():
            log.warning("发现重复基因名，添加后缀去重")
            gene_names = gene_names.to_series().pipe(
                lambda s: s.groupby(s).cumcount().astype(str).radd(
                    s.where(~s.duplicated(keep=False), s + '_')
                )
            )
            gene_names = gene_names.values

        metadata = pd.read_csv(CFG.barcodes_file, index_col=0)
        if CFG.meta_columns:
            rename_map = {}
            for target_col, source_col in CFG.meta_columns.items():
                if source_col in metadata.columns:
                    rename_map[source_col] = target_col
            if rename_map:
                metadata.rename(columns=rename_map, inplace=True)

        adata = sc.AnnData(X=mtx, obs=metadata, var=pd.DataFrame(index=gene_names))
        log.info("加载完成: %d 细胞 × %d 基因", adata.n_obs, adata.n_vars)

    elif CFG.data_format == "h5ad":
        log.info("从 h5ad 加载: %s", CFG.input_h5ad)
        adata = sc.read(CFG.input_h5ad)
        log.info("加载完成: %d 细胞 × %d 基因", adata.n_obs, adata.n_vars)

    else:
        log.error("未知 data_format: %s", CFG.data_format)
        sys.exit(1)

    # ── 保存 ──
    log.info("保存到 %s...", CFG.raw_h5ad)
    from utils import safe_write
    safe_write(adata, CFG.raw_h5ad)
    log.info("Step 00 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
