#!/usr/bin/env python3
"""
Step 06: 探索性分析
======================
  1. 细胞组成随发育阶段的变化 (堆叠柱状图)
  2. UMAP 上的 QC 指标检查
  3. 已知标记基因的 UMAP 表达
  4. 聚类大小统计

输入: 04_clustered.h5ad (需要 Stage 05 运行后以获得 cell_type 注释)
输出: CSV 表格 + PNG 图片 (不修改 h5ad)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_plot
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_composition(adata, group_col, stage_col, stage_order, fig_dir, log):
    """绘制细胞类型随发育阶段的组成变化堆积图"""
    if group_col not in adata.obs or stage_col not in adata.obs:
        log.warning("缺少 %s 或 %s，跳过组成图", group_col, stage_col)
        return
    ct_counts = adata.obs.groupby([stage_col, group_col]).size().reset_index(name='count')
    ct_pivot = ct_counts.pivot_table(
        index=stage_col, columns=group_col, values='count', fill_value=0
    )
    avail_stages = [s for s in stage_order if s in ct_pivot.index]
    if not avail_stages:
        avail_stages = list(ct_pivot.index)
    ct_pivot = ct_pivot.reindex(avail_stages)
    ct_pivot = ct_pivot.div(ct_pivot.sum(axis=1), axis=0)

    n_types = ct_pivot.shape[1]
    colors = plt.cm.tab20(np.linspace(0, 1, min(n_types, 20)))
    if n_types > 20:
        colors = np.tile(colors, int(np.ceil(n_types / 20)))[:n_types]

    fig, ax = plt.subplots(figsize=(max(10, len(avail_stages) * 1.5), 6))
    ct_pivot.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.8)
    ax.set_xlabel('Developmental stage')
    ax.set_ylabel('Fraction of cells')
    ax.set_title(f'Cluster composition by stage ({group_col})')
    ax.legend(title=group_col, bbox_to_anchor=(1.02, 1),
              loc='upper left', fontsize=8, title_fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f'composition_by_stage_{group_col}.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  组成图已保存: composition_by_stage_%s.png", group_col)

    # 导出 CSV
    ct_pivot.to_csv(os.path.join(CFG.table_dir, f'composition_by_stage_{group_col}.csv'))
    log.info("  组成表已导出")

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("06_exploratory", os.path.join(CFG.log_dir, "06_exploratory.log"))
    log.info("Step 06: 探索性分析")

    adata = sc.read(CFG.cluster_h5ad)
    log.info("加载: %s — %d 细胞", CFG.cluster_h5ad, adata.n_obs)

    fig_dir = os.path.join(CFG.figure_dir, '06_exploratory')
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir
    sc.settings.autoshow = False

    # 1. 细胞组成
    group_by = ['cell_type', 'cell_type_sub', 'leiden']
    for g in group_by:
        if g in adata.obs:
            plot_composition(adata, g, 'stage' if 'stage' in adata.obs else 'sample',
                             CFG.stage_order, fig_dir, log)

    # 2. UMAP: QC 指标
    qc_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    qc_metrics = [m for m in qc_metrics if m in adata.obs]
    if qc_metrics:
        safe_plot(sc.pl.umap, adata, color=qc_metrics, show=False,
                  save='_06_qc_umap.pdf', vmax='p99', ncols=3)

    # 3. UMAP: 标记基因
    all_markers = []
    for genes in CFG.marker_dict.values():
        all_markers.extend([g for g in genes if g in adata.raw.var_names][:2])
    all_markers = list(dict.fromkeys(all_markers))  # 去重保序
    if all_markers:
        n_markers = len(all_markers)
        batch_size = 12
        for batch_start in range(0, n_markers, batch_size):
            batch = all_markers[batch_start:batch_start + batch_size]
            safe_plot(sc.pl.umap, adata, color=batch, use_raw=True,
                      show=False, save=f'_06_markers_{batch_start}.pdf',
                      vmax='p99', ncols=4)

    # 4. 标记基因 dotplot
    if all_markers:
        safe_plot(sc.pl.dotplot, adata, var_names=all_markers,
                  groupby='cell_type', show=False, save='_06_marker_dotplot.pdf')

    # 5. 聚类大小统计
    for group_col in ['cell_type', 'leiden']:
        if group_col not in adata.obs:
            continue
        sizes = adata.obs[group_col].value_counts().sort_index()
        log.info("  %s 大小分布:", group_col)
        for label, cnt in sizes.items():
            log.info("    %s: %d cells (%.1f%%)", label, cnt, 100 * cnt / adata.n_obs)
        sizes.to_csv(os.path.join(CFG.table_dir, f'{group_col}_sizes.csv'),
                     header=['n_cells'])

    log.info("Step 06 完成, 耗时 %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
