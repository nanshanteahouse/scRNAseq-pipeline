#!/usr/bin/env python3
"""
Step 06: Interactive subclustering of a selected cell type
===========================================================
  Extract a user-selected cell type, re-run PCA + neighbors + UMAP + Leiden,
  and optionally use AI to re-annotate subclusters.

  Designed to be called in a loop (once per cell type) for iterative
  refinement of subpopulations within major cell types.

Input:  05_annotated.h5ad  (from Step 05)
Output: 05_sub_{cell_type}.h5ad  (per cell type, in CFG.h5ad_dir)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt


def main():
    t0 = time.time()

    # ── Argument parsing ──────────────────────────────────────────────
    parser = argparse.ArgumentParser(
        description="Subcluster a selected cell type from annotated data."
    )
    parser.add_argument("--config", default="../config.py",
                        help="Path to config.py")
    parser.add_argument("--cell-type", required=True,
                        help="Cell type to extract and subcluster "
                             "(e.g. 'T cell', 'Microglia')")
    args = parser.parse_args()

    # ── Config & logger ───────────────────────────────────────────────
    CFG = resolve_config(args.config)
    log = setup_logger("06_subcluster",
                       os.path.join(CFG.log_dir, "06_subcluster.log"))
    log.info("Step 06: Interactive subclustering")
    log.info("Cell type: %s", args.cell_type)

    # ── (a) Load annotated data ───────────────────────────────────────
    input_path = os.path.join(CFG.h5ad_dir, "05_annotated.h5ad")
    if not os.path.exists(input_path):
        log.error("Input file not found: %s", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d cells, %d genes",
             input_path, adata.n_obs, adata.n_vars)

    if 'cell_type' not in adata.obs:
        log.error("'cell_type' column not found in adata.obs. "
                  "Run Step 05 (annotate) first.")
        sys.exit(1)

    # ── (b) Filter to selected cell type ──────────────────────────────
    mask = adata.obs['cell_type'] == args.cell_type
    n_cells = mask.sum()
    sub = adata[mask].copy()
    log.info("Filtered: %d cells for cell type '%s'", n_cells, args.cell_type)

    # ── (c) Minimum cell check ────────────────────────────────────────
    min_cells = 50
    if n_cells < min_cells:
        log.warning("Too few cells (%d < %d). Skipping subclustering.",
                    n_cells, min_cells)
        # Still save the subset (without subcluster fields)
        safe_cell_type = args.cell_type.replace(" ", "_").replace("/", "_")
        output_path = os.path.join(
            CFG.h5ad_dir, f"05_sub_{safe_cell_type}.h5ad")
        safe_write(sub, output_path)
        log.info("Saved subset (no subclustering performed): %s", output_path)
        return

    # ── (d) PCA on subset ─────────────────────────────────────────────
    n_comps = min(50, sub.n_obs - 2)
    log.info("Running PCA (n_comps=%d)...", n_comps)
    sc.pp.pca(sub, n_comps=n_comps, svd_solver='arpack')

    # ── (e) Neighbors ─────────────────────────────────────────────────
    log.info("Computing neighbor graph (n_pcs=30)...")
    sc.pp.neighbors(sub, n_pcs=30, random_state=CFG.random_seed,
                    n_jobs=getattr(CFG, 'n_jobs', 1))

    # ── (f) UMAP ──────────────────────────────────────────────────────
    log.info("Computing UMAP...")
    sc.tl.umap(sub, random_state=CFG.random_seed)

    # ── (g) Multi-resolution Leiden ───────────────────────────────────
    log.info("Leiden subclustering, resolutions: %s", CFG.leiden_resolutions)
    for res in CFG.leiden_resolutions:
        key = f'sub_leiden_r{res}'
        sc.tl.leiden(sub, resolution=res, key_added=key,
                     random_state=CFG.random_seed)
        n_cl = sub.obs[key].nunique()
        log.info("  r=%.1f → %d subclusters", res, n_cl)

    # ── (h) Set best-resolution leiden ─────────────────────────────────
    best_key = f'sub_leiden_r{CFG.best_resolution}'
    if best_key in sub.obs:
        sub.obs['leiden'] = sub.obs[best_key].copy()
        log.info("  Best resolution: sub_leiden_r%.1f → leiden set",
                 CFG.best_resolution)
    else:
        avail = [k for k in sub.obs if k.startswith('sub_leiden_')]
        if avail:
            sub.obs['leiden'] = sub.obs[avail[-1]].copy()
            log.info("  Fallback to %s for 'leiden'", avail[-1])
        else:
            log.warning("No Leiden results available — skipping cluster label.")
            sub.obs['leiden'] = "0"

    # ── (i) Save UMAP plots ───────────────────────────────────────────
    sc.settings.figdir = CFG.figure_dir
    sc.settings.autoshow = False
    safe_cell_type = args.cell_type.replace(" ", "_").replace("/", "_")

    # Leiden at best resolution
    safe_plot(sc.pl.umap, sub, color='leiden', show=False,
              save=f'_06_sub_{safe_cell_type}_leiden.pdf',
              legend_loc='on data', title=f'{args.cell_type} — leiden')

    # Multi-resolution comparison
    res_keys = [f'sub_leiden_r{r}' for r in CFG.leiden_resolutions
                if f'sub_leiden_r{r}' in sub.obs]
    i = -1
    if res_keys:
        n_res = len(res_keys)
        n_cols = min(3, n_res)
        n_rows = int(np.ceil(n_res / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols,
                                 figsize=(6 * n_cols, 5 * n_rows))
        axes = axes.ravel() if n_res > 1 else [axes]
        for i, key in enumerate(res_keys):
            sc.pl.umap(sub, color=key, ax=axes[i], show=False,
                       legend_loc='on data', legend_fontsize=6,
                       title=f'{safe_cell_type} — {key}')
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')
        fig.tight_layout()
        fig.savefig(os.path.join(CFG.figure_dir,
                                 f'umap_sub_{safe_cell_type}_resolutions.pdf'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        log.info("Multi-resolution UMAP saved for %s", safe_cell_type)

    # ── (j) AI-based subcluster annotation ────────────────────────────
    if CFG.ai.enabled and CFG.ai.ai_subcluster:
        try:
            from ai_prompts import build_annotation_prompt
            from ai_caller import ai_query

            log.info("Running AI subcluster re-annotation...")

            # build_annotation_prompt runs rank_genes_groups internally
            # and returns (system_prompt, user_prompt)
            sys_prompt, user_prompt = build_annotation_prompt(
                sub, tissue=args.cell_type, species="unknown",
            )

            result = ai_query(sys_prompt, user_prompt, CFG.ai)
            log.info("AI response received (%d chars)", len(result))

            # ── Parse JSON from AI response ──
            # Strip potential markdown code fences
            cleaned = result.strip()
            if cleaned.startswith("```"):
                lines = cleaned.split("\n")
                start = 0
                for i, line in enumerate(lines):
                    if line.strip().startswith("```"):
                        start = i + 1
                        break
                end = len(lines)
                for i in range(len(lines) - 1, start - 1, -1):
                    if lines[i].strip().startswith("```"):
                        end = i
                        break
                cleaned = "\n".join(lines[start:end]).strip()

            parsed = json.loads(cleaned)

            # Build per-cluster labels with cell_type + subtype
            ai_labels = {}
            for cluster_id, info in parsed.items():
                cell_type = info.get("cell_type", "Unknown")
                subtype = info.get("subtype", "N/A")
                if subtype and subtype.upper() != "N/A":
                    label = f"{cell_type}_{subtype}"
                else:
                    label = cell_type
                # Sanitize for categorical use
                label = label.replace(" ", "_").replace("/", "_")
                ai_labels[cluster_id] = label

            # Map to sub.obs as categorical
            sub.obs['sub_ai_label'] = (
                sub.obs['leiden'].astype(str).map(ai_labels)
            ).astype('category')
            n_ai_types = sub.obs['sub_ai_label'].nunique()
            log.info("AI annotation: %d subcluster types identified",
                     n_ai_types)

            # Log per-cluster AI mapping
            for cluster_id in sorted(sub.obs['leiden'].unique(),
                                     key=lambda x: int(x)):
                label = ai_labels.get(str(cluster_id), "Unmapped")
                count = (sub.obs['leiden'] == cluster_id).sum()
                log.info("  Subcluster %s → %s (%d cells)",
                         cluster_id, label, count)

            # Save AI-annotated UMAP
            safe_plot(sc.pl.umap, sub, color='sub_ai_label', show=False,
                      save=f'_06_sub_{safe_cell_type}_ai.pdf',
                      legend_loc='on data',
                      title=f'{args.cell_type} — AI subcluster')

        except Exception as e:
            log.warning("AI subcluster annotation failed: %s", e)
            log.info("Falling back to numeric subcluster labels.")
            sub.obs['sub_ai_label'] = (
                "Subcluster_" + sub.obs['leiden'].astype(str)
            ).astype('category')
    else:
        log.info("AI subcluster annotation disabled "
                 "(CFG.ai.enabled=%s, CFG.ai.ai_subcluster=%s)",
                 CFG.ai.enabled, CFG.ai.ai_subcluster)

    # ── (k) Save subcluster result ────────────────────────────────────
    safe_filename = f"05_sub_{safe_cell_type}.h5ad"
    output_path = os.path.join(CFG.h5ad_dir, safe_filename)
    safe_write(sub, output_path)
    log.info("Saved: %s", output_path)

    # ── (l) Summary ───────────────────────────────────────────────────
    n_clusters = sub.obs['leiden'].nunique()
    log.info("=" * 50)
    log.info("Subcluster Summary")
    log.info("  Cell type:       %s", args.cell_type)
    log.info("  Cells:           %d", n_cells)
    log.info("  Subclusters:     %d", n_clusters)
    log.info("  Resolution:      %.1f", CFG.best_resolution)
    log.info("  Per-cluster counts:")
    for cluster_id in sorted(sub.obs['leiden'].unique(),
                             key=lambda x: int(x)):
        count = (sub.obs['leiden'] == cluster_id).sum()
        pct = 100.0 * count / n_cells
        log.info("    Cluster %s: %d cells (%.1f%%)",
                 cluster_id, count, pct)
    log.info("=" * 50)
    log.info("Step 06 done, %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
