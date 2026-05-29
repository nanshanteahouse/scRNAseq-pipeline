#!/usr/bin/env python3
"""
Step 04: 邻居图 + UMAP + 多分辨率 Leiden 聚类
==================================================
  - 在 Harmony 校正后的 PCA 上建图
  - 多分辨率聚类 (GSE138002 策略: 6 个分辨率)
  - 保存所有分辨率结果用于后续选择

输入: 03_harmony.h5ad
输出: 04_clustered.h5ad (含邻居图, UMAP, 多分辨率 Leiden 标签)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("04_cluster", os.path.join(CFG.log_dir, "04_cluster_umap.log"))
    log.info("Step 04: 邻居图 + UMAP + Leiden 聚类")

    adata = sc.read(CFG.harmony_h5ad)
    log.info("加载: %s", adata.shape)

    use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'

    # ── 邻居图 ──
    log.info("计算邻居图 (n_neighbors=%d, use_rep=%s)...",
             CFG.n_neighbors, use_rep)
    sc.pp.neighbors(
        adata, n_neighbors=CFG.n_neighbors,
        n_pcs=CFG.n_pcs_use, use_rep=use_rep,
        random_state=CFG.random_seed,
    )

    # ── UMAP ──
    log.info("计算 UMAP...")
    sc.tl.umap(adata, min_dist=0.3, spread=1.0, random_state=CFG.random_seed)

    # 预可视化
    for color_by in ['sample', 'stage']:
        if color_by in adata.obs:
            safe_plot(sc.pl.umap, adata, color=color_by, show=False,
                      save=f'_04_{color_by}.pdf',
                      title=f'UMAP colored by {color_by}')

    # ── 多分辨率 Leiden 聚类 ──
    log.info("Leiden 聚类, 分辨率: %s", CFG.leiden_resolutions)
    for res in CFG.leiden_resolutions:
        key = f'leiden_r{res}'
        sc.tl.leiden(adata, resolution=res, key_added=key,
                      random_state=CFG.random_seed,
                      flavor=CFG.leiden_flavor)
        n_clusters = adata.obs[key].nunique()
        log.info("  r=%.1f → %d clusters", res, n_clusters)

    # 设置主分辨率
    best_key = f'leiden_r{CFG.best_resolution}'
    if best_key in adata.obs:
        adata.obs['leiden'] = adata.obs[best_key].copy()
        log.info("  主分辨率: leiden_r%.1f", CFG.best_resolution)
    else:
        # 回退到最后一个可用分辨率
        avail = [k for k in adata.obs if k.startswith('leiden_')]
        if avail:
            adata.obs['leiden'] = adata.obs[avail[-1]].copy()
            log.info("  回退到 %s", avail[-1])

    # 多分辨率 UMAP 对比图
    res_keys = [f'leiden_r{r}' for r in CFG.leiden_resolutions
                if f'leiden_r{r}' in adata.obs]
    n_res = len(res_keys)
    if n_res > 0:
        n_cols = min(3, n_res)
        n_rows = int(np.ceil(n_res / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
        axes = axes.ravel() if n_res > 1 else [axes]
        for i, key in enumerate(res_keys):
            sc.pl.umap(adata, color=key, ax=axes[i], show=False,
                       legend_loc='on data', legend_fontsize=6,
                       title=key)
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')
        fig.tight_layout()
        fig.savefig(os.path.join(CFG.figure_dir, 'umap_leiden_resolutions.pdf'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        log.info("  多分辨率 UMAP 图已保存")

    safe_write(adata, CFG.cluster_h5ad)
    log.info("Step 04 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
