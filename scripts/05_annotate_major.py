#!/usr/bin/env python3
"""
Step 05: 细胞类型自动注释 (Major Lineage — AI + Score_genes 双模式)
=====================================================================
  双模式注释策略:
    1. AI 模式 (首选): 基于 marker 基因 + LLM 智能判断细胞类型
    2. Score_genes 模式 (回退): 基于已知 marker 基因打分
  输出主要细胞类型 (major lineage) 及亚型/状态/置信度信息。

输入: 04_clustered.h5ad
输出: 05_annotated.h5ad (新增 cell_type, cell_subtype, cell_state, annot_confidence, ... 列)
"""
import sys, os, time, argparse, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import pandas as pd
import numpy as np
import logging

log: logging.Logger


# ═══════════════════════════════════════════════════════════════════════
#  旧有注释函数 (Score_genes 模式)
# ═══════════════════════════════════════════════════════════════════════

def run_annotation(adata, marker_dict, logger):
    """基于 marker 基因打分的细胞类型注释 (来自 05_annotate.py 原有逻辑)。"""
    if not marker_dict:
        logger.warning("未配置 marker_dict，跳过注释。")
        adata.obs['cell_type'] = adata.obs['leiden'].astype(str)
        return

    cell_types = list(marker_dict.keys())
    logger.info("评分注释: %d 种待选类型", len(cell_types))

    for ct in cell_types:
        genes = marker_dict[ct]
        genes_present = [g for g in genes if g in adata.raw.var_names]
        if not genes_present:
            logger.warning("  %s: 无 marker 基因在数据中", ct)
            adata.obs[f'score_{ct}'] = 0.0
            continue
        sc.tl.score_genes(adata, gene_list=genes_present,
                          score_name=f'score_{ct}', random_state=42)

    # 每个聚类取最高分的类型
    score_cols = [f'score_{ct}' for ct in cell_types]
    groupby_kw = {'observed': True} if hasattr(pd.Categorical, 'observed') else {}
    cluster_scores = adata.obs.groupby('leiden', **groupby_kw)[score_cols].mean()
    best_match = cluster_scores.idxmax(axis=1)
    best_ct = best_match.str.replace('score_', '')

    cluster_to_ct = dict(zip(best_ct.index, best_ct.values))
    adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_ct).astype('category')

    logger.info("聚类 → 细胞类型映射:")
    for label in sorted(adata.obs['leiden'].unique()):
        ct = cluster_to_ct[label]
        max_score = cluster_scores.loc[label, f'score_{ct}']
        logger.info("  聚类 %s → %s (score=%.3f)", label, ct, max_score)

    # 置信度: 最高分与次高分之差
    if len(cell_types) >= 2:
        sorted_scores = cluster_scores.apply(
            lambda row: row.sort_values(ascending=False).values, axis=1, result_type='expand'
        )
        confidence = sorted_scores.iloc[:, 0] - sorted_scores.iloc[:, 1]
        adata.obs['annotation_confidence'] = adata.obs['leiden'].map(confidence).astype(float).values
        low_conf = (adata.obs['annotation_confidence'] < 0.02).sum()
        if low_conf > 0:
            logger.info("  低置信度细胞 (<0.02): %d (%.1f%%)",
                        low_conf, 100 * low_conf / adata.n_obs)

    logger.info("注释完成: %d 种细胞类型", adata.obs['cell_type'].nunique())


def run_subclustering(adata, subcluster_types, resolution, min_cells, logger):
    """基于 parent cell_type 的子聚类 (来自 05_annotate.py 原有逻辑)。"""
    if not subcluster_types:
        logger.info("未配置子聚类类型，跳过。")
        adata.obs['cell_type_sub'] = adata.obs['cell_type'].astype(str)
        return

    logger.info("子聚类: %s (resolution=%.1f)...", subcluster_types, resolution)
    adata.obs['cell_type_sub'] = adata.obs['cell_type'].astype(str)

    for parent_type in subcluster_types:
        mask = adata.obs['cell_type'] == parent_type
        n_cells = mask.sum()
        if n_cells < min_cells:
            logger.info("  %s: 细胞太少 (%d < %d), 跳过", parent_type, n_cells, min_cells)
            continue

        logger.info("  子聚类 %s (%d cells)...", parent_type, n_cells)
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
    logger.info("子聚类完成: %d 种亚型", n_sub)


# ═══════════════════════════════════════════════════════════════════════
#  AI 注释函数
# ═══════════════════════════════════════════════════════════════════════

def ai_annotate(adata, CFG, logger):
    """
    基于 LLM 的 AI 注释主流程。

    步骤:
      1. rank_genes_groups → 获取各聚类 marker 基因
      2. 保存 marker 基因 CSV
      3. 构建提示词 → 调用 LLM
      4. 解析 JSON 响应 → 映射注释到 adata.obs
      5. 生成 UMAP 可视化 & 导出注释表格

    返回:
        annotations dict (解析成功) 或 None (失败，触发回退)
    """
    # ── a. 计算 marker 基因 ───────────────────────────────────────────
    logger.info("计算 marker 基因 (Wilcoxon rank-sum)...")
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

    # ── b. 保存 marker 基因 CSV ───────────────────────────────────────
    marker_rows = []
    for cl in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
        df = sc.get.rank_genes_groups_df(adata, group=str(cl))
        df['cluster'] = cl
        marker_rows.append(df)
    marker_df = pd.concat(marker_rows, ignore_index=True)
    marker_csv = os.path.join(CFG.table_dir, 'marker_genes_per_group.csv')
    marker_df.to_csv(marker_csv, index=False)
    logger.info("Marker 基因已保存: %s", marker_csv)

    # ── c. 获取组织 & 物种 ────────────────────────────────────────────
    tissue = CFG.tissue
    species = CFG.species
    logger.info("注释上下文: tissue=%s, species=%s", tissue, species)

    # ── d. 构建提示词 ─────────────────────────────────────────────────
    from ai_prompts import build_annotation_prompt
    sys_prompt, user_prompt = build_annotation_prompt(adata, tissue, species, precomputed_rank=True)

    # ── e. 调用 LLM ───────────────────────────────────────────────────
    from ai_caller import ai_query
    logger.info("向 LLM 请求细胞类型注释 (model=%s)...", CFG.ai.model)
    try:
        response = ai_query(sys_prompt, user_prompt, cfg=CFG.ai)
    except Exception as exc:
        logger.warning("LLM 查询失败: %s — 回退到 score_genes 方法", exc)
        return None

    # ── f. 解析 JSON ──────────────────────────────────────────────────
    try:
        annotations = json.loads(response)
    except json.JSONDecodeError:
        logger.warning("LLM 响应不是有效 JSON — 回退到 score_genes 方法")
        logger.warning("原始响应 (前 500 字符): %s", response[:500])
        return None

    # 验证每聚类注释结构
    required_keys = {'cell_type', 'state', 'subtype', 'confidence', 'reasoning'}
    for cid, ann in annotations.items():
        if not isinstance(ann, dict):
            logger.warning("聚类 %s 注释不是 dict — 回退", cid)
            return None
        missing = required_keys - ann.keys()
        if missing:
            logger.warning("聚类 %s 缺少字段 %s — 回退", cid, missing)
            return None

    logger.info("LLM 注释成功: %d 个聚类已解析", len(annotations))

    # ── g/h. 映射注释到 adata.obs ─────────────────────────────────────
    leiden_str = adata.obs['leiden'].astype(str)
    adata.obs['cell_type'] = leiden_str.map(
        {k: v['cell_type'] for k, v in annotations.items()}
    ).astype('category')
    adata.obs['cell_state'] = leiden_str.map(
        {k: v['state'] for k, v in annotations.items()}
    )
    adata.obs['cell_subtype'] = leiden_str.map(
        {k: v['subtype'] for k, v in annotations.items()}
    )
    adata.obs['annot_confidence'] = leiden_str.map(
        {k: v['confidence'] for k, v in annotations.items()}
    )
    adata.obs['annot_reasoning'] = leiden_str.map(
        {k: v['reasoning'] for k, v in annotations.items()}
    )

    # ── i. 保存注释 CSV ───────────────────────────────────────────────
    ann_records = []
    for cid in sorted(annotations.keys(), key=lambda x: int(x)):
        ann = annotations[cid]
        ann_records.append({
            'cluster': cid,
            'cell_type': ann['cell_type'],
            'state': ann['state'],
            'subtype': ann['subtype'],
            'confidence': ann['confidence'],
            'reasoning': ann['reasoning'],
        })
    ann_df = pd.DataFrame(ann_records)
    ann_csv = os.path.join(CFG.table_dir, 'cell_type_annotations.csv')
    ann_df.to_csv(ann_csv, index=False)
    logger.info("注释表已保存: %s", ann_csv)

    # 日志输出映射
    logger.info("聚类 → 细胞类型映射 (AI):")
    for rec in ann_records:
        logger.info("  聚类 %s → %s (state=%s, subtype=%s, conf=%s)",
                    rec['cluster'], rec['cell_type'],
                    rec['state'], rec['subtype'], rec['confidence'])

    # ── j. UMAP 可视化 ────────────────────────────────────────────────
    sc.settings.figdir = CFG.figure_dir
    sc.settings.autoshow = False

    # annot_label = cell_type (+ state 如果不为 N/A)
    adata.obs['annot_label'] = adata.obs['cell_type'].astype(str)
    state_not_na = adata.obs['cell_state'] != 'N/A'
    adata.obs.loc[state_not_na, 'annot_label'] = (
        adata.obs.loc[state_not_na, 'cell_type'].astype(str)
        + ' (' + adata.obs.loc[state_not_na, 'cell_state'].astype(str) + ')'
    )

    safe_plot(sc.pl.umap, adata, color='cell_type', show=False,
              save='_05_celltype_ai.pdf', legend_loc='on data')
    safe_plot(sc.pl.umap, adata, color='annot_label', show=False,
              save='_05_annot_label_ai.pdf', legend_loc='on data')

    # ── k. 细胞元数据导出 ─────────────────────────────────────────────
    meta_df = pd.DataFrame({
        'barcode': adata.obs_names,
        'UMAP_1': adata.obsm['X_umap'][:, 0],
        'UMAP_2': adata.obsm['X_umap'][:, 1],
        'cell_type': adata.obs['cell_type'].values,
        'cell_state': adata.obs['cell_state'].values,
        'cell_subtype': adata.obs['cell_subtype'].values,
        'annot_confidence': adata.obs['annot_confidence'].values,
    })
    meta_csv = os.path.join(CFG.table_dir, 'cell_metadata.csv')
    meta_df.to_csv(meta_csv, index=False)
    logger.info("细胞元数据已导出: %s", meta_csv)

    return annotations


# ═══════════════════════════════════════════════════════════════════════
#  Score_genes 模式 (回退)
# ═══════════════════════════════════════════════════════════════════════

def score_genes_mode(adata, CFG, logger):
    """Score_genes 回退模式 — 复用旧有 run_annotation + run_subclustering。"""
    logger.info("Score_genes 模式 — 基于标记基因打分注释")

    run_annotation(adata, CFG.marker_dict, logger)
    run_subclustering(adata, CFG.subcluster_types,
                      CFG.subcluster_resolution, CFG.min_cells_subcluster, logger)

    # 统一列名: cell_type_sub → cell_subtype
    if 'cell_type_sub' in adata.obs:
        adata.obs['cell_subtype'] = adata.obs['cell_type_sub'].astype(str)

    # annot_label (这里仅为 cell_type，无 state 信息)
    adata.obs['annot_label'] = adata.obs['cell_type'].astype(str)

    # 置信度重命名
    if 'annotation_confidence' in adata.obs:
        adata.obs['annot_confidence'] = adata.obs['annotation_confidence']

    # 可视化
    sc.settings.figdir = CFG.figure_dir
    sc.settings.autoshow = False
    safe_plot(sc.pl.umap, adata, color='cell_type', show=False,
              save='_05_celltype.pdf', legend_loc='on data')
    safe_plot(sc.pl.umap, adata, color='annot_label', show=False,
              save='_05_annot_label.pdf', legend_loc='on data')
    if 'annotation_confidence' in adata.obs:
        safe_plot(sc.pl.umap, adata, color='annotation_confidence', show=False,
                  save='_05_confidence.pdf', cmap='viridis')

    # 导出细胞元数据
    meta_cols = ['barcode']
    if 'X_umap' in adata.obsm:
        meta_df = pd.DataFrame({
            'barcode': adata.obs_names,
            'UMAP_1': adata.obsm['X_umap'][:, 0],
            'UMAP_2': adata.obsm['X_umap'][:, 1],
        })
    else:
        meta_df = pd.DataFrame({'barcode': adata.obs_names})
    for col in ['cell_type', 'cell_subtype', 'cell_type_sub', 'annotation_confidence']:
        if col in adata.obs:
            meta_df[col] = adata.obs[col].values
    meta_csv = os.path.join(CFG.table_dir, 'cell_metadata.csv')
    meta_df.to_csv(meta_csv, index=False)
    logger.info("细胞元数据已导出: %s", meta_csv)


# ═══════════════════════════════════════════════════════════════════════
#  主入口
# ═══════════════════════════════════════════════════════════════════════

def main():
    global log
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("05_annotate_major",
                        os.path.join(CFG.log_dir, "05_annotate_major.log"))
    log.info("Step 05: 细胞类型注释 (Major Lineage)")

    adata = sc.read(CFG.cluster_h5ad)
    log.info("加载: %s — %d 细胞, %d 聚类",
             CFG.cluster_h5ad, adata.n_obs, adata.obs['leiden'].nunique())

    # ── 判断 AI 模式是否可用 ──────────────────────────────────────────
    ai_enabled = getattr(CFG.ai, 'enabled', False)
    ai_annot_on = getattr(CFG.ai, 'ai_annotation', False)

    if ai_enabled and ai_annot_on:
        log.info("AI 模式启用 — 使用 LLM 进行智能注释")
        ann_result = ai_annotate(adata, CFG, log)
        if ann_result is not None:
            # AI 注释成功 → 保存并退出
            safe_write(adata, CFG.annotated_h5ad)
            log.info("Step 05 (AI mode) 完成, 耗时 %.1fs", time.time() - t0)
            return
        # AI 注释失败 → 记录并回退
        log.warning("AI 注释失败，回退到 Score_genes 模式")

    # ── Score_genes 模式 (AI 未启用 或 回退) ──────────────────────────
    score_genes_mode(adata, CFG, log)
    safe_write(adata, CFG.annotated_h5ad)
    log.info("Step 05 (score_genes mode) 完成, 耗时 %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
