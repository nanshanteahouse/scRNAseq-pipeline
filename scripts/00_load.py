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
import scipy.sparse as sp

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
        # Legacy 2-column genes.tsv.gz → 3-column features.tsv.gz
        genes_path = os.path.join(CFG.mtx_dir, CFG.mtx_prefix + 'genes.tsv.gz')
        features_path = os.path.join(CFG.mtx_dir, CFG.mtx_prefix + 'features.tsv.gz')
        if not os.path.exists(features_path) and os.path.exists(genes_path):
            log.info("检测到旧版 2-column genes.tsv.gz — 转换为 features.tsv.gz...")
            import gzip
            with gzip.open(genes_path, 'rt') as f_in:
                with gzip.open(features_path, 'wt') as f_out:
                    for line in f_in:
                        f_out.write(line.rstrip('\n') + '\tGene Expression\n')
            log.info("  features.tsv.gz 已创建")

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
        matrix_ext = os.path.splitext(CFG.matrix_file)[1].lower()
        if matrix_ext in ('.csv', '.gz'):
            # True CSV format: gene × cell, first column = gene names
            log.info("从 CSV 加载: %s", CFG.matrix_file)
            df = pd.read_csv(CFG.matrix_file, index_col=0)
            log.info("CSV 形状: %s", df.shape)
            # Transpose to AnnData convention: cells × genes
            adata = sc.AnnData(X=df.values.T.astype(np.float32))
            adata.var_names = df.index.astype(str)
            adata.obs_names = df.columns.astype(str)
            # Load metadata if barcodes/features files provided
            if CFG.barcodes_file and os.path.exists(CFG.barcodes_file):
                metadata = pd.read_csv(CFG.barcodes_file, index_col=0)
                adata.obs = adata.obs.join(metadata, how='left')
            if CFG.features_file and os.path.exists(CFG.features_file):
                genes = pd.read_csv(CFG.features_file)
                if len(genes) == adata.n_vars:
                    adata.var = genes
        else:
            # Original MTX path (mmread)
            log.info("从 MTX 矩阵加载: %s", CFG.matrix_file)
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
        backed = getattr(CFG, 'backed', None) or None
        adata = sc.read(CFG.input_h5ad, backed=backed) if backed else sc.read(CFG.input_h5ad)
        log.info("加载完成: %d 细胞 × %d 基因", adata.n_obs, adata.n_vars)

    else:
        log.error("未知 data_format: %s", CFG.data_format)
        sys.exit(1)

    # ── 统一稀疏格式: CSR (行优先) ──
    if getattr(CFG, 'force_csr', True) and sp.issparse(adata.X):
        if not sp.isspmatrix_csr(adata.X):
            adata.X = adata.X.tocsr()
            log.info("X 格式已转换为 CSR")

    # ── 可选 float32 精度 ──
    if getattr(CFG, 'use_float32', False):
        adata.X = adata.X.astype('float32', copy=False) if sp.issparse(adata.X) else adata.X
        log.info("X 精度已转换为 float32")

    # ── 保存 ──
    log.info("保存到 %s...", CFG.raw_h5ad)
    from utils import safe_write
    safe_write(adata, CFG.raw_h5ad)
    log.info("Step 00 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
