#!/usr/bin/env python3
"""
Step 04: 邻居图 + UMAP + 多参数网格 Leiden 聚类
==================================================
  - 在 Harmony 校正后的 PCA 上建图
  - 多参数网格扫描 (n_neighbors × resolution)
  - 保存所有组合结果用于交互比较

输入: 03_integrated.h5ad
输出: 04_grid_results.h5ad (含所有参数组合的邻居图、UMAP、Leiden 标签)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("04_cluster", os.path.join(CFG.log_dir, "04_cluster_umap.log"))
    log.info("Step 04: 邻居图 + UMAP + 多参数网格 Leiden 聚类")

    # ── 输入 ──
    input_path = CFG.integrated_h5ad
    log.info("加载: %s", input_path)
    adata = sc.read(input_path)
    log.info("  shape: %s", adata.shape)

    use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'
    log.info("use_rep: %s", use_rep)

    # ── 参数网格 ──
    n_neighbors_grid = getattr(CFG, 'param_grid_n_neighbors', [15, 20, 30])
    resolutions_grid = getattr(CFG, 'param_grid_resolutions', [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
    log.info("参数网格: n_neighbors=%s, resolutions=%s", n_neighbors_grid, resolutions_grid)

    results_summary = []

    for n in n_neighbors_grid:
        # 邻居图
        log.info("计算邻居图 (n_neighbors=%d, use_rep=%s)...", n, use_rep)
        try:
            sc.pp.neighbors(
                adata, n_neighbors=n,
                n_pcs=CFG.n_pcs_use, use_rep=use_rep,
                random_state=CFG.random_seed,
            )
        except Exception as e:
            log.error("邻居图计算失败 (n_neighbors=%d): %s", n, e)
            continue

        for res in resolutions_grid:
            umap_key = f'umap_{n}_{res}'
            leiden_key = f'leiden_{n}_{res}'

            log.info("  UMAP + Leiden (n=%d, r=%.1f)...", n, res)
            try:
                # UMAP
                sc.tl.umap(adata, min_dist=0.3, spread=1.0,
                           random_state=CFG.random_seed)
                # 保留当前 UMAP 坐标
                adata.obsm[umap_key] = adata.obsm['X_umap'].copy()

                # Leiden 聚类
                sc.tl.leiden(adata, resolution=res, key_added=leiden_key,
                             random_state=CFG.random_seed,
                             flavor=CFG.leiden_flavor)
            except Exception as e:
                log.error("  UMAP/Leiden 失败 (n=%d, r=%.1f): %s", n, res, e)
                continue

            n_clusters = adata.obs[leiden_key].nunique()
            log.info("    → %d clusters", n_clusters)

            # Silhouette score (在 PCA 空间上计算，大数据集抽样)
            sil_score = None
            try:
                pca_key = use_rep
                if adata.n_obs > 10000:
                    rng = np.random.RandomState(CFG.random_seed)
                    idx = rng.choice(adata.n_obs, 10000, replace=False)
                    sil_score = silhouette_score(
                        adata.obsm[pca_key][idx, :CFG.n_pcs_use],
                        adata.obs[leiden_key].values[idx],
                    )
                else:
                    sil_score = silhouette_score(
                        adata.obsm[pca_key][:, :CFG.n_pcs_use],
                        adata.obs[leiden_key].values,
                    )
                log.info("    silhouette_score=%.4f", sil_score)
            except Exception as e:
                log.warning("    silhouette_score 计算失败: %s", e)

            results_summary.append({
                'n_neighbors': n,
                'resolution': res,
                'n_clusters': n_clusters,
                'silhouette_score': sil_score,
            })

            # 单参数组合 UMAP 图
            try:
                safe_plot(sc.pl.umap, adata, color=leiden_key, show=False,
                          title=f'UMAP (n_neighbors={n}, resolution={res})')
                plt.savefig(
                    os.path.join(CFG.figure_dir,
                                 f'umap_grid_n{n}_r{res}.png'),
                    dpi=150, bbox_inches='tight')
                plt.close()
                log.info("    图已保存: umap_grid_n%d_r%.1f.png", n, res)
            except Exception as e:
                log.warning("    单参数 UMAP 图保存失败: %s", e)

    # ── 汇总 CSV ──
    df_summary = pd.DataFrame(results_summary)
    csv_path = os.path.join(CFG.table_dir, 'param_grid_summary.csv')
    try:
        df_summary.to_csv(csv_path, index=False)
        log.info("参数网格汇总已保存: %s", csv_path)
        log.info("\n%s", df_summary.to_string())
    except Exception as e:
        log.warning("汇总 CSV 保存失败: %s", e)

    # ── 网格汇总图: 所有参数组合对比 ──
    n_n = len(n_neighbors_grid)
    n_r = len(resolutions_grid)
    try:
        fig, axes = plt.subplots(n_n, n_r,
                                 figsize=(5 * n_r + 2, 4 * n_n + 1),
                                 squeeze=False)
        for i, n in enumerate(n_neighbors_grid):
            for j, res in enumerate(resolutions_grid):
                ax = axes[i, j]
                umap_key = f'umap_{n}_{res}'
                leiden_key = f'leiden_{n}_{res}'
                if umap_key in adata.obsm and leiden_key in adata.obs:
                    saved_umap = adata.obsm['X_umap'].copy()
                    try:
                        adata.obsm['X_umap'] = adata.obsm[umap_key].copy()
                        sc.pl.umap(adata, color=leiden_key, ax=ax,
                                   show=False, legend_loc='on data',
                                   legend_fontsize=5,
                                   title=f'n={n}, r={res}')
                    except Exception as e_sub:
                        log.warning("  子图失败 (n=%d, r=%.1f): %s",
                                    n, res, e_sub)
                        ax.text(0.5, 0.5, 'Error', ha='center',
                                va='center', transform=ax.transAxes)
                    finally:
                        adata.obsm['X_umap'] = saved_umap
                else:
                    ax.text(0.5, 0.5, 'N/A', ha='center', va='center',
                            transform=ax.transAxes, fontsize=12)
                    ax.set_title(f'n={n}, r={res}')
        fig.tight_layout()
        fig.savefig(os.path.join(CFG.figure_dir, 'umap_param_grid_summary.png'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        log.info("参数网格汇总图已保存")
    except Exception as e:
        log.warning("网格汇总图生成失败: %s", e)

    # ── 按 n_neighbors 分组的多分辨率对比图 ──
    for n in n_neighbors_grid:
        res_keys = [f'leiden_{n}_{r}' for r in resolutions_grid
                    if f'leiden_{n}_{r}' in adata.obs]
        n_res = len(res_keys)
        if n_res > 0:
            try:
                n_cols = min(3, n_res)
                n_rows = int(np.ceil(n_res / n_cols))
                fig, axes = plt.subplots(n_rows, n_cols,
                                         figsize=(6 * n_cols, 5 * n_rows))
                axes = axes.ravel() if n_res > 1 else [axes]
                for i, key in enumerate(res_keys):
                    sc.pl.umap(adata, color=key, ax=axes[i], show=False,
                               legend_loc='on data', legend_fontsize=6,
                               title=key)
                for j in range(len(res_keys), len(axes)):
                    axes[j].axis('off')
                fig.tight_layout()
                fig.savefig(
                    os.path.join(CFG.figure_dir,
                                 f'umap_leiden_n{n}_all_resolutions.pdf'),
                    dpi=150, bbox_inches='tight')
                plt.close(fig)
                log.info("  多分辨率 UMAP 图 (n=%d) 已保存", n)
            except Exception as e:
                log.warning("  多分辨率对比图 (n=%d) 失败: %s", n, e)

    # ── 保存临时 h5ad (非最终 checkpoint) ──
    temp_path = os.path.join(CFG.h5ad_dir, "04_grid_results.h5ad")
    try:
        safe_write(adata, temp_path)
        log.info("临时 h5ad 已保存: %s", temp_path)
    except Exception as e:
        log.error("临时 h5ad 保存失败: %s", e)

    log.info("Step 04 完成, 耗时 %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
