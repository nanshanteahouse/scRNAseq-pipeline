#!/usr/bin/env python3
"""
Step 03: PCA 降维 + Harmony 批次校正
========================================
  - 先用 n_pcs_full (100) 算 PCA 供 elbow 参考
  - Harmony 在校正后的嵌入存入 X_pca_harmony
  - 不覆盖原始 X_pca (保留校正前后对比能力)

输入: 02_normalized.h5ad
输出: 03_harmony.h5ad (含 X_pca, X_pca_harmony)
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
    log = setup_logger("03_pca_harmony", os.path.join(CFG.log_dir, "03_pca_harmony.log"))
    log.info("Step 03: PCA + Harmony 批次校正")

    adata = sc.read(CFG.norm_h5ad)
    log.info("加载: %s", adata.shape)

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

    safe_write(adata, CFG.harmony_h5ad)
    log.info("Step 03 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
