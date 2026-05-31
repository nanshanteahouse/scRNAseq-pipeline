#!/usr/bin/env python3
"""
Step 03: 归一化 + HVG 选择 + PCA + Harmony 批次校正（整合版）
===============================================================
整合原 Step 02 + Step 03 为单一步骤，并加入 regress_out。

关键顺序 (继承 GSE169109 最佳实践):
  1. raw counts → regress_out (total_counts, pct_counts_mt)
  2. raw counts → 找 HVG (seurat_v3, batch-aware)
  3. normalize_total → log1p
  4. adata.raw = adata  (保留全基因)
  5. X 只保留 HVG 用于下游降维
  6. PCA (n_pcs_full, elbow 图)
  7. Harmony 批次校正

输入: 02_qc.h5ad
输出: 03_integrated.h5ad (X = log1p(normalized) on HVGs, .raw = 全基因,
                          obsm: X_pca, X_pca_harmony)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import matplotlib.pyplot as plt

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("03_integrate", os.path.join(CFG.log_dir, "03_integrate.log"))
    log.info("Step 03: 归一化 + HVG 选择 + PCA + Harmony（整合版）")

    # ── 读取 ──
    adata = sc.read(CFG.qc_h5ad)
    log.info("加载: %d 细胞 × %d 基因", adata.n_obs, adata.n_vars)

    # ── Regress out 技术变异 ──
    # 在 raw counts 上回归 total_counts 和 pct_counts_mt，移除技术噪音
    # 某些数据集可能没有 pct_counts_mt，用 try/except 兜底
    try:
        log.info("回归技术协变量: total_counts, pct_counts_mt ...")
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        log.info("  regress_out 完成")
    except Exception as e:
        log.warning("regress_out 失败（跳过）: %s", e)

    # ── HVG ──
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

    # ── 归一化 ──
    log.info("归一化 (target_sum=%.0f) + log1p...", CFG.normalize_target_sum)
    sc.pp.normalize_total(adata, target_sum=CFG.normalize_target_sum)
    sc.pp.log1p(adata)

    # ── 保留全基因到 .raw ──
    adata.raw = adata
    log.info(".raw 已保存 (全基因: %d vars)", adata.raw.n_vars)

    # ── X 只保留 HVG ──
    adata = adata[:, adata.var['highly_variable']].copy()
    log.info("X 缩小到 HVGs: %s", adata.shape)

    # ── PCA ──
    log.info("PCA (%d components)...", CFG.n_pcs_full)
    sc.pp.pca(adata, n_comps=CFG.n_pcs_full,
              svd_solver='arpack', random_state=CFG.random_seed)
    var_ratio = adata.uns['pca']['variance_ratio']
    log.info("  top-5 方差比: %.4f", var_ratio[:5].sum())
    log.info("  前 50 PC 累积方差比: %.4f", var_ratio[:50].sum())

    # PCA elbow 图
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(1, CFG.n_pcs_full + 1), var_ratio, 'o-', ms=3)
    ax.axvline(CFG.n_pcs_use, color='red', linestyle='--', alpha=0.5,
               label=f'n_pcs_use={CFG.n_pcs_use}')
    ax.set_xlabel('PC'); ax.set_ylabel('Variance ratio')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(CFG.figure_dir, 'pca_elbow.png'), dpi=150)
    plt.close(fig)
    log.info("  PCA elbow 图已保存")

    # ── Harmony ──
    if CFG.use_harmony:
        from harmony import harmonize
        batch_key = CFG.harmony_batch_key
        if batch_key not in adata.obs:
            log.warning("Harmony batch_key '%s' 不在 obs 中，跳过校正", batch_key)
        else:
            log.info("Harmony 校正 (batch_key=%s)...", batch_key)
            Z = harmonize(
                adata.obsm['X_pca'][:, :CFG.n_pcs_use],
                adata.obs,
                batch_key=batch_key,
                random_state=CFG.random_seed,
                max_iter_harmony=CFG.harmony_max_iter,
            )
            adata.obsm['X_pca_harmony'] = Z
            log.info("  Harmony 完成, 输出形状: %s", Z.shape)
            # 对比图
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            sc.pl.embedding(adata, basis='X_pca', color=batch_key,
                            ax=axes[0], show=False, title='PCA (before Harmony)')
            sc.pl.embedding(adata, basis='X_pca_harmony', color=batch_key,
                            ax=axes[1], show=False, title='Harmony-corrected')
            fig.tight_layout()
            fig.savefig(os.path.join(CFG.figure_dir, 'harmony_comparison.png'), dpi=150)
            plt.close(fig)
            log.info("  Harmony 对比图已保存")
    else:
        log.info("Harmony 已禁用，使用原始 PCA。")
        adata.obsm['X_pca_harmony'] = adata.obsm['X_pca'].copy()

    # ── 保存 ──
    out_path = os.path.join(CFG.h5ad_dir, "03_integrated.h5ad")
    safe_write(adata, out_path)
    log.info("Step 03 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
