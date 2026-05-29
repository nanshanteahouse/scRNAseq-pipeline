#!/usr/bin/env python3
"""
Step 05: 细胞类型自动注释 + 子聚类
========================================
  混合策略 (GSE169109 + GSE138002):
    1. 对每个聚类计算已知标记基因的得分
    2. 每个聚类分配到得分最高的细胞类型
    3. 对主要类型做子聚类 (如 NP/Neuron)
    4. 输出标注结果和置信度

输入: 04_clustered.h5ad
输出: 同文件 (新增 cell_type, cell_type_sub, score_* 列)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_plot
import scanpy as sc
import pandas as pd
import numpy as np

def run_annotation(adata, marker_dict, log):
    if not marker_dict:
        log.warning("未配置 marker_dict，跳过注释。")
        adata.obs['cell_type'] = adata.obs['leiden'].astype(str)
        return

    cell_types = list(marker_dict.keys())
    log.info("评分注释: %d 种待选类型", len(cell_types))

    for ct in cell_types:
        genes = marker_dict[ct]
        genes_present = [g for g in genes if g in adata.raw.var_names]
        if not genes_present:
            log.warning("  %s: 无 marker 基因在数据中", ct)
            adata.obs[f'score_{ct}'] = 0.0
            continue
        sc.tl.score_genes(adata, gene_list=genes_present,
                          score_name=f'score_{ct}', random_state=42)

    # 每个聚类取最高分的类型
    score_cols = [f'score_{ct}' for ct in cell_types]
    # 使用 observed=True 避免 FutureWarning (category 类型)
    groupby_kw = {'observed': True} if hasattr(pd.Categorical, 'observed') else {}
    cluster_scores = adata.obs.groupby('leiden', **groupby_kw)[score_cols].mean()
    best_match = cluster_scores.idxmax(axis=1)
    best_ct = best_match.str.replace('score_', '')

    cluster_to_ct = dict(zip(best_ct.index, best_ct.values))
    adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_ct).astype('category')

    log.info("聚类 → 细胞类型映射:")
    for label in sorted(adata.obs['leiden'].unique()):
        ct = cluster_to_ct[label]
        max_score = cluster_scores.loc[label, f'score_{ct}']
        log.info("  聚类 %s → %s (score=%.3f)", label, ct, max_score)

    # 置信度: 最高分与次高分之差
    if len(cell_types) >= 2:
        sorted_scores = cluster_scores.apply(
            lambda row: row.sort_values(ascending=False).values, axis=1, result_type='expand'
        )
        confidence = sorted_scores.iloc[:, 0] - sorted_scores.iloc[:, 1]
        adata.obs['annotation_confidence'] = adata.obs['leiden'].map(confidence).values
        low_conf = (adata.obs['annotation_confidence'] < 0.02).sum()
        if low_conf > 0:
            log.info("  低置信度细胞 (<0.02): %d (%.1f%%)",
                     low_conf, 100 * low_conf / adata.n_obs)

    log.info("注释完成: %d 种细胞类型", adata.obs['cell_type'].nunique())

def run_subclustering(adata, subcluster_types, resolution, min_cells, log):
    if not subcluster_types:
        log.info("未配置子聚类类型，跳过。")
        adata.obs['cell_type_sub'] = adata.obs['cell_type'].astype(str)
        return

    log.info("子聚类: %s (resolution=%.1f)...", subcluster_types, resolution)
    adata.obs['cell_type_sub'] = adata.obs['cell_type'].astype(str)

    for parent_type in subcluster_types:
        mask = adata.obs['cell_type'] == parent_type
        n_cells = mask.sum()
        if n_cells < min_cells:
            log.info("  %s: 细胞太少 (%d < %d), 跳过", parent_type, n_cells, min_cells)
            continue

        log.info("  子聚类 %s (%d cells)...", parent_type, n_cells)
        sub = adata[mask].copy()
        sc.pp.neighbors(sub, n_pcs=50, use_rep='X_pca_harmony',
                        random_state=42)
        sc.tl.leiden(sub, resolution=resolution, key_added='subcluster',
                     random_state=42)
        labels = np.array(sub.obs['cell_type'].astype(str)
                          + '_' + sub.obs['subcluster'].astype(str))
        adata.obs.loc[mask, 'cell_type_sub'] = labels.tolist()

    adata.obs['cell_type_sub'] = adata.obs['cell_type_sub'].astype('category')
    n_sub = adata.obs['cell_type_sub'].nunique()
    log.info("子聚类完成: %d 种亚型", n_sub)

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("05_annotate", os.path.join(CFG.log_dir, "05_annotate.log"))
    log.info("Step 05: 细胞类型注释")

    adata = sc.read(CFG.cluster_h5ad)
    log.info("加载: %s — %d 细胞, %d 聚类",
             CFG.cluster_h5ad, adata.n_obs, adata.obs['leiden'].nunique())

    run_annotation(adata, CFG.marker_dict, log)
    run_subclustering(adata, CFG.subcluster_types,
                      CFG.subcluster_resolution, CFG.min_cells_subcluster, log)

    # 可视化
    sc.settings.figdir = CFG.figure_dir
    sc.settings.autoshow = False
    safe_plot(sc.pl.umap, adata, color='cell_type', show=False,
              save='_05_celltype.pdf', legend_loc='on data')
    safe_plot(sc.pl.umap, adata, color='cell_type_sub', show=False,
              save='_05_subtype.pdf', legend_loc='on data')
    safe_plot(sc.pl.umap, adata, color='annotation_confidence', show=False,
              save='_05_confidence.pdf', cmap='viridis')

    # 更新保存
    from utils import safe_write
    safe_write(adata, CFG.cluster_h5ad)
    log.info("Step 05 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
