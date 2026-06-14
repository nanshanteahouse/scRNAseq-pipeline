#!/usr/bin/env python3
"""
Step 09: GO/KEGG 富集分析
=============================
输入: Step 07 输出的 marker_genes_per_group.csv
      （每类细胞 vs 其他所有细胞的 Wilcoxon 标记基因）

方法 (通过 GSEApy + Enrichr API):
  ORA: 取每类上调基因 top N → 过表达分析
  Pre-ranked GSEA: 使用全部基因的 score 排序 → 无需 cutoff

输出:
  tables/enrichment/
    ora_{gene_set}_{cluster}.csv         — ORA 结果表
    prerank_{gene_set}_{cluster}.csv     — GSEA 结果表
    ora_{gene_set}_summary.csv           — 汇总（所有聚类合并）
  figures/enrichment/
    ora_{gene_set}_bubble.pdf            — 气泡图
    ora_{gene_set}_dotplot.pdf           — 点图

依赖: pip install gseapy (需要 Rust 编译器)
      curl https://sh.rustup.rs -sSf | sh  # 如需要
"""
import json
import sys, os, time, argparse, warnings
from typing import Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)


def read_marker_csv(table_dir: str, log) -> pd.DataFrame:
    """读取 Step 07 产出的标记基因 CSV"""
    path = os.path.join(table_dir, "marker_genes_per_group.csv")
    if not os.path.exists(path):
        log.error("未找到标记基因文件: %s", path)
        log.error("请先运行 Step 07 (07_markers_de.py) 生成该文件。")
        sys.exit(1)
    df = pd.read_csv(path)
    log.info("加载标记基因: %d 行, %d 个分组",
             len(df), df['group'].nunique())
    return df


def _ora_one_group(
    grp,
    grp_df: pd.DataFrame,
    gene_set: str,
    CFG,
    log,
):
    """Run ORA for a single group via Enrichr API (used by ThreadPoolExecutor)."""
    import gseapy as gp

    grp_up = (
        grp_df[grp_df['logfoldchanges'] > 0]
        .nsmallest(CFG.enrichment_n_top_genes, 'pvals_adj')
        ['names']
        .str.upper()
        .dropna()
        .unique()
        .tolist()
    )
    if len(grp_up) < CFG.enrichment_min_size:
        log.info("  %s: 上调基因不足 (%d < %d), 跳过 ORA",
                 grp, len(grp_up), CFG.enrichment_min_size)
        return (grp, None)

    try:
        enr = gp.enrichr(
            gene_list=grp_up,
            gene_sets=gene_set,
            organism=CFG.enrichment_organism,
            outdir=None,
            no_plot=True,
            verbose=False,
        )
    except Exception as e:
        log.warning("  %s ORA 失败 (%s): %s", grp, gene_set, e)
        return (grp, None)

    res = enr.results
    if res is None or len(res) == 0:
        return (grp, None)
    res['cluster'] = grp
    res['n_genes_input'] = len(grp_up)
    n_sig = (res['Adjusted P-value'] < CFG.enrichment_pval_cutoff).sum()
    log.info("  %s: %d/%d 显著通路 (ORA, %s)",
             grp, n_sig, len(res), gene_set.split('_')[0])
    return (grp, res)


def run_ora(
    marker_df: pd.DataFrame,
    gene_set: str,
    CFG,
    log,
) -> pd.DataFrame:
    """
    ORA (Over-Representation Analysis) via GSEApy Enrichr.

    对每个分组的 top N 上调基因（按 pvals_adj 排序），
    查询 Enrichr API 计算通路富集。
    """
    import gseapy as gp

    groups = marker_df['group'].unique()
    max_workers = min(5, getattr(CFG, 'n_jobs', 4))
    results = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_grp = {
            executor.submit(
                _ora_one_group, grp,
                marker_df[marker_df['group'] == grp],
                gene_set, CFG, log,
            ): grp
            for grp in groups
        }
        for future in as_completed(future_to_grp):
            grp, res = future.result()
            results[grp] = res

    all_rows = [results[grp] for grp in groups if results.get(grp) is not None]

    if not all_rows:
        log.warning("  ORA 无结果 (gene_set=%s)", gene_set)
        return pd.DataFrame()

    combined = pd.concat(all_rows, ignore_index=True)
    return combined


def _prerank_one_group(
    grp,
    grp_df: pd.DataFrame,
    gene_set: str,
    CFG,
    log,
):
    """Run pre-ranked GSEA for a single group (used by ThreadPoolExecutor)."""
    import gseapy as gp

    grp_df = grp_df.dropna(subset=['scores', 'names'])
    if len(grp_df) < CFG.enrichment_min_size:
        log.info("  %s: 基因不足 (%d < %d), 跳过 GSEA",
                 grp, len(grp_df), CFG.enrichment_min_size)
        return (grp, None)

    rnk = grp_df.set_index('names')['scores'].drop_duplicates()
    rnk.index = rnk.index.str.upper()

    try:
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_set,
            min_size=CFG.enrichment_min_size,
            max_size=CFG.enrichment_max_size,
            permutation_num=CFG.enrichment_permutations,
            threads=1,
            outdir=None,
            seed=CFG.random_seed,
            verbose=False,
            no_plot=True,
        )
    except Exception as e:
        log.warning("  %s GSEA 失败 (%s): %s", grp, gene_set, e)
        return (grp, None)

    res = pre_res.res2d
    if res is None or len(res) == 0:
        return (grp, None)
    res['cluster'] = grp
    res['n_genes_input'] = len(rnk)
    n_sig = (res['FDR q-val'] < CFG.enrichment_pval_cutoff).sum()
    log.info("  %s: %d/%d 显著通路 (GSEA, %s)",
             grp, n_sig, len(res), gene_set.split('_')[0])
    return (grp, res)


def run_prerank(
    marker_df: pd.DataFrame,
    gene_set: str,
    CFG,
    log,
) -> pd.DataFrame:
    """
    Pre-ranked GSEA via GSEApy.

    对每个分组使用全部基因的 scores 作为排序指标，
    无需 cutoff，捕获微弱的协同变化。
    """
    import gseapy as gp

    groups = marker_df['group'].unique()
    max_workers = min(5, getattr(CFG, 'n_jobs', 4))
    results = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_grp = {
            executor.submit(
                _prerank_one_group, grp,
                marker_df[marker_df['group'] == grp],
                gene_set, CFG, log,
            ): grp
            for grp in groups
        }
        for future in as_completed(future_to_grp):
            grp, res = future.result()
            results[grp] = res

    all_rows = [results[grp] for grp in groups if results.get(grp) is not None]

    if not all_rows:
        log.warning("  GSEA 无结果 (gene_set=%s)", gene_set)
        return pd.DataFrame()

    combined = pd.concat(all_rows, ignore_index=True)
    return combined


def save_results(
    ora_results: dict,
    prerank_results: dict,
    CFG,
    log,
) -> None:
    """保存富集结果 CSV + 绘图"""
    table_dir = os.path.join(CFG.table_dir, "enrichment")
    fig_dir = os.path.join(CFG.figure_dir, "enrichment")
    os.makedirs(table_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    # ── 保存 ORA CSV ──
    for gs_name, df in ora_results.items():
        if df.empty:
            continue
        path = os.path.join(table_dir, f"ora_{gs_name}_summary.csv")
        df.to_csv(path, index=False)
        log.info("  ORA 结果导出: %s (%d 行)", path, len(df))

    # ── 保存 GSEA CSV ──
    for gs_name, df in prerank_results.items():
        if df.empty:
            continue
        path = os.path.join(table_dir, f"prerank_{gs_name}_summary.csv")
        df.to_csv(path, index=False)
        log.info("  GSEA 结果导出: %s (%d 行)", path, len(df))

    # ── 气泡图（每组 top 20 通路） ──
    try:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        for gs_name, df in ora_results.items():
            if df.empty:
                continue
            plot_ora_bubble(df, gs_name, fig_dir, CFG, log)

        for gs_name, df in prerank_results.items():
            if df.empty:
                continue
            plot_prerank_bubble(df, gs_name, fig_dir, CFG, log)
    except Exception as e:
        log.warning("绘图失败: %s", e)


def plot_ora_bubble(
    df: pd.DataFrame,
    gs_name: str,
    fig_dir: str,
    CFG,
    log,
) -> None:
    """ORA 气泡图: x=cluster, y=Term, size=Overlap, color=Adjusted P-value"""
    import matplotlib.pyplot as plt

    sig = df[df['Adjusted P-value'] < CFG.enrichment_pval_cutoff].copy()
    if sig.empty:
        sig = df.head(5)
    top_per_cluster = (
        sig.sort_values('Adjusted P-value')
        .groupby('cluster')
        .head(CFG.enrichment_n_top_genes // max(1, sig['cluster'].nunique()))
    )
    if len(top_per_cluster) < 3:
        log.info("  跳过气泡图 (%s): 显著通路不足", gs_name)
        return

    top_per_cluster['-log10_padj'] = -np.log10(
        top_per_cluster['Adjusted P-value'].clip(lower=1e-300)
    )
    # 简短 Term 名称
    top_per_cluster['Term_short'] = top_per_cluster['Term'].str.replace(
        r'\s*\(GO:\d+\)$', '', regex=True
    ).str[:60]

    fig, ax = plt.subplots(figsize=(
        max(8, 0.5 * top_per_cluster['cluster'].nunique()),
        max(6, 0.3 * top_per_cluster['Term_short'].nunique()),
    ))
    overlap_numeric = (
        top_per_cluster['Overlap']
        .astype(str).str.split('/').str[0].astype(float)
    )
    sc = ax.scatter(
        top_per_cluster['cluster'],
        top_per_cluster['Term_short'],
        s=overlap_numeric * 30,
        c=top_per_cluster['-log10_padj'],
        cmap='YlOrRd',
        edgecolors='grey', linewidths=0.5,
    )
    plt.colorbar(sc, ax=ax, label='-log10(Adjusted P-value)')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('')
    ax.set_title(f'Enrichment: {gs_name}')
    fig.tight_layout()
    path = os.path.join(fig_dir, f"ora_{gs_name}_bubble.pdf")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  ORA 气泡图: %s", path)


def plot_prerank_bubble(
    df: pd.DataFrame,
    gs_name: str,
    fig_dir: str,
    CFG,
    log,
) -> None:
    """Pre-ranked GSEA 气泡图: color=NES, size=-log10(FDR)"""
    import matplotlib.pyplot as plt

    sig = df[df['FDR q-val'] < CFG.enrichment_pval_cutoff].copy()
    if sig.empty:
        sig = df.head(10)
    top_per_cluster = (
        sig.sort_values('FDR q-val')
        .groupby('cluster')
        .head(CFG.enrichment_n_top_genes // max(1, sig['cluster'].nunique()))
    )
    if len(top_per_cluster) < 3:
        return

    top_per_cluster['Term_short'] = top_per_cluster['Term'].str[:60]
    top_per_cluster['-log10_fdr'] = -np.log10(
        top_per_cluster['FDR q-val'].clip(lower=1e-300)
    )
    # NES 颜色: 红色=上调, 蓝色=下调
    vmax = max(abs(top_per_cluster['NES'].min()),
               abs(top_per_cluster['NES'].max()))
    norm = plt.Normalize(-vmax, vmax)

    fig, ax = plt.subplots(figsize=(
        max(8, 0.5 * top_per_cluster['cluster'].nunique()),
        max(6, 0.3 * top_per_cluster['Term_short'].nunique()),
    ))
    sc = ax.scatter(
        top_per_cluster['cluster'],
        top_per_cluster['Term_short'],
        s=top_per_cluster['-log10_fdr'] * 20,
        c=top_per_cluster['NES'],
        cmap='RdBu_r', norm=norm,
        edgecolors='grey', linewidths=0.5,
    )
    plt.colorbar(sc, ax=ax, label='NES')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('')
    ax.set_title(f'GSEA: {gs_name}')
    fig.tight_layout()
    path = os.path.join(fig_dir, f"prerank_{gs_name}_bubble.pdf")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  GSEA 气泡图: %s", path)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("09_enrichment",
                        os.path.join(CFG.log_dir, "09_enrichment.log"))
    log.info("Step 09: 富集分析 (GO/KEGG)")

    if not CFG.run_enrichment:
        log.info("富集分析已禁用 (run_enrichment=False)")
        return

    marker_df = read_marker_csv(CFG.table_dir, log)
    log.info("基因集库: %s", CFG.enrichment_gene_sets)
    log.info("方法: %s", CFG.enrichment_method)

    # ── 按基因集库循环 ──
    ora_results = {}
    prerank_results = {}

    for gs in CFG.enrichment_gene_sets:
        gs_name = gs.replace(' ', '_').replace('/', '_')

        if CFG.enrichment_method in ('ora', 'both'):
            log.info("[ORA] 基因集: %s", gs)
            ora_df = run_ora(marker_df, gs, CFG, log)
            if ora_df is not None and not ora_df.empty:
                ora_results[gs_name] = ora_df

        if CFG.enrichment_method in ('prerank', 'both'):
            log.info("[GSEA] 基因集: %s", gs)
            prerank_df = run_prerank(marker_df, gs, CFG, log)
            if prerank_df is not None and not prerank_df.empty:
                prerank_results[gs_name] = prerank_df

    total_ora = sum(len(df) for df in ora_results.values())
    total_gsea = sum(len(df) for df in prerank_results.values())
    log.info("富集结果汇总: ORA %d 行, GSEA %d 行",
             total_ora, total_gsea)

    save_results(ora_results, prerank_results, CFG, log)

    # ── AI Biological Interpretation (optional) ──
    if CFG.ai.enabled and CFG.ai.ai_interpretation:
        log.info("AI: Generating biological interpretation...")
        try:
            summary_data = []
            for gs_name, df in ora_results.items():
                if df.empty:
                    continue
                sig = df[df['Adjusted P-value'] < CFG.enrichment_pval_cutoff]
                for cluster in sig['cluster'].unique():
                    cluster_sig = sig[sig['cluster'] == cluster].head(5)
                    summary_data.append({
                        "gene_set": gs_name,
                        "cluster": str(cluster),
                        "top_terms": cluster_sig[['Term', 'Adjusted P-value']].to_dict('records')
                    })

            if summary_data:
                system_prompt = "You are an expert computational biologist interpreting scRNA-seq enrichment results."
                user_prompt = f"Enrichment results summary:\n{json.dumps(summary_data, indent=2)}\n\nProvide biological interpretation: key pathways, cross-cell-type patterns, and testable hypotheses."

                from scripts.ai_caller import ai_query
                interpretation = ai_query(system_prompt, user_prompt, cfg=CFG.ai)

                interp_path = os.path.join(CFG.table_dir, "enrichment", "ai_interpretation.txt")
                os.makedirs(os.path.dirname(interp_path), exist_ok=True)
                with open(interp_path, "w") as f:
                    f.write(interpretation)
                log.info("AI interpretation saved to %s", interp_path)

                summary_lines = [f"Biological Interpretation — {'scRNA-seq enrichment'}"]
                summary_lines.append("=" * 60)
                summary_lines.append(interpretation[:2000])
                summary_path = os.path.join(CFG.table_dir, "enrichment", "ai_interpretation_summary.txt")
                with open(summary_path, "w") as f:
                    f.write("\n".join(summary_lines))
        except Exception as e:
            log.warning("AI interpretation skipped: %s", e)

    log.info("Step 09 完成, 耗时 %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
