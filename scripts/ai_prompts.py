#!/usr/bin/env python3
"""
ai_prompts.py — AI 注释与解读的提示词模板
=============================================

集中存放所有 LLM 提示词模板和构建函数，确保提示词一致、可复用、易维护。

设计原则:
  - 提示词与调用逻辑分离（提示词在此模块，调用在 ai_caller.py）
  - ANNOTATION_SYSTEM_PROMPT 为独立常量，可直接导入使用
  - ANNOTATION_USER_PROMPT_TEMPLATE 使用 format() 占位符
  - build_annotation_prompt() 自动运行 marker 基因检测并组装完整提示词

用法:
    from ai_prompts import ANNOTATION_SYSTEM_PROMPT, build_annotation_prompt
    sys_prompt, user_prompt = build_annotation_prompt(adata, "hypothalamus", "human")
    result = ai_query(sys_prompt, user_prompt, cfg)
"""

import json
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import scanpy as sc


# ═══════════════════════════════════════════════════════════════════════
#  聚类注释提示词
# ═══════════════════════════════════════════════════════════════════════

ANNOTATION_SYSTEM_PROMPT = """You are an expert single-cell RNA-seq biologist with deep knowledge of cell type identification across tissues and species.

For each cluster ID provided in the user message, return a JSON object mapping cluster IDs to annotations with the following fields:
  - cell_type  : the broad cell type (e.g. "T cell", "Macrophage", "Oligodendrocyte", "Excitatory neuron")
  - state      : activation or functional state (e.g. "resting", "activated", "proliferating", "N/A")
  - subtype    : the most specific subtype (e.g. "CD8+ cytotoxic T cell", "M1 macrophage", "SST+ interneuron", "N/A")
  - confidence : one of "high", "medium", or "low" — based on how specific and well-established the markers are
  - reasoning  : a single sentence citing the key marker genes that support your annotation

Return ONLY a valid JSON object. No explanation, no markdown formatting, no code fences.
Include ALL cluster IDs in the response.

Required format:
{"0":{"cell_type":"T cell","state":"activated","subtype":"CD8+ cytotoxic T cell","confidence":"high","reasoning":"High CD8A, GZMB, PRF1 expression indicates cytotoxic T cells"},"1":{"cell_type":"...","state":"...","subtype":"...","confidence":"...","reasoning":"..."}}"""


ANNOTATION_USER_PROMPT_TEMPLATE = """Tissue: {tissue}
Species: {species}

Marker genes per cluster (top 20 by Wilcoxon score):
{cluster_markers_json}

Return ONLY a valid JSON object mapping each cluster ID to its annotation. Include ALL cluster IDs."""


def build_annotation_prompt(adata, tissue: str, species: str,
                            precomputed_rank: bool = False,
                            extra_context: str = "",
                            compact: bool = False):
    """
    构建聚类注释的完整提示词对。

    可自动运行或跳过 Wilcoxon rank-sum 检验。当调用者已经执行过
    rank_genes_groups 时，传入 precomputed_rank=True 避免重复计算，
    直接使用 adata.uns['rank_genes_groups'] 中的已有结果。

    参数:
        adata:  已聚类（leiden 列）的 AnnData 对象
        tissue: 组织名称（如 "hypothalamus", "retina"）
        species: 物种名称（如 "human", "mouse"）
        precomputed_rank: 若为 True，跳过 rank_genes_groups 计算，
                          使用 adata 中已有的结果（默认 False）
        extra_context: 额外上下文信息（如发育阶段列表）追加到用户提示词尾部
        compact: 若为 True，每聚类仅展示 top 10 而非 top 20 marker 基因，
                 用于减少 tokens（默认 False）

    返回:
        (system_prompt, user_prompt) 二元组，可直接传入 ai_query()
    """
    # ── 计算 marker 基因（如尚未计算）────────────────────────────────
    if not precomputed_rank:
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    # ── 提取每聚类 marker 基因 ───────────────────────────────────────
    n_top = 10 if compact else 20
    clusters = sorted(adata.obs["leiden"].unique(),
                      key=lambda x: int(x))
    cluster_markers: dict = {}
    for cl in clusters:
        df = sc.get.rank_genes_groups_df(adata, group=str(cl))
        top_genes = df.head(n_top)["names"].tolist()
        cluster_markers[cl] = top_genes

    # ── 组装提示词 ────────────────────────────────────────────────────
    user_prompt = ANNOTATION_USER_PROMPT_TEMPLATE.format(
        tissue=tissue,
        species=species,
        cluster_markers_json=json.dumps(cluster_markers, indent=2),
    )
    if extra_context:
        user_prompt += f"\n\n{extra_context}"

    return ANNOTATION_SYSTEM_PROMPT, user_prompt


# ═══════════════════════════════════════════════════════════════════════
#  后续步骤的提示词桩（TODO: 在对应步骤实现时完善）
# ═══════════════════════════════════════════════════════════════════════

# PARAM_SUGGEST_PROMPT = """..."""
# 用途: 根据数据特征建议 QC 参数阈值
# 输入: QC 统计摘要（n_genes, pct_mito, complexity 分布）
# 输出: JSON { "min_genes": ..., "max_genes": ..., "max_pct_mito": ... }

# QC_REVIEW_PROMPT = """..."""
# 用途: 审查 QC 结果并给出质量判断
# 输入: 过滤前后的细胞/基因统计对比
# 输出: JSON { "verdict": "pass"|"warn"|"fail", "issues": [...], "recommendation": "..." }

# DEG_DESIGN_PROMPT = """..."""
# 用途: 建议差异表达分析的对比设计
# 输入: 样本元数据、感兴趣的变量
# 输出: JSON { "contrasts": [{"group_a": ..., "group_b": ..., "label": "..."}] }

# INTERPRETATION_PROMPT = """..."""
# 用途: 解读差异表达或富集分析结果
# 输入: DEG 表格或富集结果
# 输出: 一段自然语言总结
