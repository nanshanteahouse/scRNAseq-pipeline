# scRNA-seq 管线输出文件参考

> 输出文件说明详细参考 — 从 README 迁移至此

本文档完整记录了管线运行后生成的各类输出文件：h5ad 检查点链、差异表达结果表、轨迹分析产出、富集分析结果，以及所有可视化图片。

---

## 目录结构

所有输出统一存放在 `results/` 目录下，按类型分三个子目录：

```
results/
├── h5ad/                              # AnnData 检查点文件
│   ├── 00_raw.h5ad                    # Step 00 — 原始数据加载
│   ├── 01_doublet.h5ad                # Step 02 — 双细胞检测后
│   ├── 02_qc.h5ad                     # Step 03 — QC 过滤后
│   ├── 03_integrated.h5ad             # Step 04 — Harmony 批次校正后
│   ├── 04_clustered.h5ad              # Step 05 — 聚类 + UMAP 后
│   ├── 05_annotated.h5ad              # Step 06 — 细胞类型注释后
│   └── 05_final.h5ad                  # Step 08 — 轨迹分析后（含 PAGA + DPT）
│
├── figures/                           # 可视化图片
│   ├── pca_elbow.png                  # Step 04 — PCA 肘部图（用于选主成分数）
│   ├── harmony_comparison.png         # Step 04 — Harmony 校正前后对比
│   ├── umap_leiden_resolutions.pdf    # Step 05 — 多分辨率 UMAP 对比图
│   ├── 05_celltype.pdf                # Step 06 — 注释后的细胞类型 UMAP
│   ├── 06_exploratory/                # Step 11 — 探索性分析图集
│   │   ├── composition_by_stage_*.png #   各阶段细胞类型组成堆叠图
│   │   ├── *_sizes.csv                #   聚类大小统计
│   │   ├── _06_qc_umap.pdf            #   UMAP 上 QC 指标分布
│   │   ├── _06_markers_*.pdf          #   已知标记基因 UMAP 表达
│   │   └── _06_marker_dotplot.pdf     #   标记基因点图
│   ├── 07_marker_heatmap.pdf          # Step 08 — 标记基因热图（每类 top5）
│   ├── 07_dotplot.pdf                 # Step 08 — 已知标记基因点图
│   ├── 08_pseudotime.pdf              # Step 09 — 扩散伪时间 UMAP
│   ├── 08_pseudotime_diffmap.pdf      # Step 09 — 扩散图上的伪时间
│   ├── 08_paga_graph.pdf              # Step 09 — PAGA 轨迹图
│   ├── 08_paga_umap.pdf               # Step 09 — UMAP 上的 PAGA 覆盖
│   ├── 08_dev_genes_heatmap.pdf       # Step 09 — 发育基因沿伪时间热图
│   ├── 08_trend_*.pdf                 # Step 09 — 单基因沿伪时间散点图
│   ├── 08_final_*.pdf                 # Step 09 — 最终 UMAP（stage/celltype/pseudotime）
│   └── enrichment/                    # Step 10 — 富集分析图集
│       ├── ora_{gene_set}_bubble.pdf  #   ORA 气泡图（size=overlap, color=P-value）
│       └── prerank_{gene_set}_bubble.pdf  # GSEA 气泡图（color=NES, size=FDR）
│
└── tables/                            # 表格输出（CSV）
    ├── marker_genes_per_group.csv      # Step 08 — Layer 1 标记基因
    ├── marker_genes_per_group_{col}.csv # 按多注释层级分别输出
    ├── pairwise_stage_de.csv           # Step 08 — Layer 2 阶段配对 DE
    ├── temporal_trend_genes.csv        # Step 08 — Layer 3 时间趋势基因
    ├── branch_deg.csv                  # Step 09 — 轨迹分支 DEG
    ├── composition_by_stage_*.csv      # Step 11 — 各阶段细胞组成表
    ├── cell_type_sizes.csv             # Step 11 — 各类型细胞数
    ├── leiden_sizes.csv                # Step 11 — 各聚类细胞数
    ├── enrichment_ora.csv              # Step 10 — ORA 汇总（兼容入口）
    ├── enrichment_gsea.csv             # Step 10 — GSEA 汇总（兼容入口）
    └── enrichment/                     # Step 10 — 详细富集结果
        ├── ora_{gene_set}_summary.csv  #   各基因集 ORA 结果表
        ├── prerank_{gene_set}_summary.csv # 各基因集 GSEA 结果表
        └── ai_interpretation.txt       #   AI 生物学解读报告
```

> **注意**：`figures/enrichment/` 和 `tables/enrichment/` 中的 `{gene_set}` 占位符对应配置中的基因集库名称，例如 `GO_Biological_Process`、`KEGG_2021_Human`、`Reactome_2022` 等。

---

## h5ad 检查点链

每个 h5ad 文件是一个独立完整的 AnnData 对象，包含该步骤的全部分析结果。这种设计使得：

- **断点恢复**：`--resume` 模式扫描检查点缺失情况，自动从第一个未完成的步骤继续
- **独立验证**：任意 h5ad 可单独加载验证，无需重跑前面步骤
- **灵活回溯**：发现某步参数不佳时，只需从该步骤重新开始

```python
import scanpy as sc
adata = sc.read("results/h5ad/03_integrated.h5ad")
adata  # 包含 raw、PCA、Harmony 校正后的嵌入
```

各检查点的核心内容：

| 文件 | 包含 | 典型用途 |
|------|------|----------|
| `00_raw.h5ad` | 原始 counts、样本/阶段元数据 | 数据完整性校验 |
| `01_doublet.h5ad` | Scrublet 双细胞分数 + 预测标签 | 双细胞率统计 |
| `02_qc.h5ad` | QC 指标、过滤掩码 | QC 阈值审查 |
| `03_integrated.h5ad` | 归一化表达量、HVG、PCA、Harmony 校正嵌入 | 降维质量检查 |
| `04_clustered.h5ad` | UMAP、Leiden 聚类标签、多分辨率结果 | 聚类结构探索 |
| `05_annotated.h5ad` | 细胞类型注释、AI 注释置信度 | 注释结果审查 |
| `05_final.h5ad` | PAGA 图、扩散图、DPT 伪时间、分支结果 | 轨迹分析审查 |

---

## 差异表达输出（DE）

管线设计了三层差异表达分析（Step 08 `07_markers_de.py`）和一层轨迹分支差异表达（Step 09 `08_trajectory.py`），从不同角度揭示基因表达变化。

### Layer 1：标记基因检测

**文件**：`tables/marker_genes_per_group.csv`

每类细胞 vs 其余所有细胞的差异表达（Wilcoxon rank-sum 检验），用于鉴定各类细胞的特征标记基因。

| 列名 | 含义 |
|------|------|
| `group` | 细胞类型（或聚类编号）名称 |
| `names` | 基因符号 |
| `scores` |  Wilcoxon 秩和统计量 |
| `logfoldchanges` | 平均表达倍数的 log2 变化 |
| `pvals` | 原始 P 值 |
| `pvals_adj` | Bonferroni 校正后 P 值 |
| `pts` | 该组中表达该基因的细胞比例 |
| `pts_rest` | 其余组中表达该基因的细胞比例 |

如果数据中有多个注释层级（`cell_type`、`cell_type_sub`、`leiden`），每个层级会单独输出一个 `marker_genes_per_group_{column}.csv` 文件。主层级同时输出为 `marker_genes_per_group.csv` 作为兼容入口。

### Layer 2：相邻阶段配对 DE

**文件**：`tables/pairwise_stage_de.csv`

对同一类细胞，比较其在相邻发育阶段之间的表达差异（t 检验）。用于追踪细胞类型在发育过程中的转录变化。

| 列名 | 含义 |
|------|------|
| `cell_type` | 细胞类型 |
| `comparison` | 比较对，格式 `后阶段_vs_前阶段` |
| `names` | 基因符号 |
| `scores` | 检验统计量 |
| `logfoldchanges` | 倍数变化 |
| `pvals_adj` | 校正后 P 值 |

只有当数据中包含 `stage` 注释且配置了 `stage_order` 时才运行。每对阶段要求每组至少 5 个细胞。

### Layer 3：发育时间趋势基因

**文件**：`tables/temporal_trend_genes.csv`

对每类细胞，计算各基因的平均表达量随发育阶段的 Spearman 相关。筛选出表达量单调上升或下降的基因（各取前 20）。

| 列名 | 含义 |
|------|------|
| `cell_type` | 细胞类型 |
| `gene` | 基因符号 |
| `spearman_r` | Spearman 相关系数（正=上升，负=下降） |
| `direction` | `up`（上升）或 `down`（下降） |

要求至少 3 个阶段且每阶段至少 5 个细胞。每个方向取 top 20 基因。

### 分支 DEG

**文件**：`tables/branch_deg.csv`

在 PAGA 轨迹中，对分支决策点的两个子谱系进行差异表达分析。用于鉴定驱动细胞命运决定的谱系特化基因。

| 列名 | 含义 |
|------|------|
| `branch` | 分支对，格式 `子谱系_vs_父谱系` |
| `names` | 基因符号 |
| `scores` | 检验统计量 |
| `logfoldchanges` | 倍数变化 |
| `pvals_adj` | 校正后 P 值 |

---

## 富集分析输出

### 方法概述

富集分析（Step 10 `09_enrichment.py`）通过 GSEApy 调用 Enrichr Web API，支持两种互补方法：

**ORA（Over-Representation Analysis，过表达分析）**
- 输入：每类细胞的上调标记基因（按校正 P 值取 top N）
- 方法：超几何检验（Fisher 精确检验的变体）
- 特点：结果直观，富集/不富集二值判断
- 输出文件：`enrichment_ora.csv` 及 `tables/enrichment/ora_*_summary.csv`

**Pre-ranked GSEA（基因集富集分析）**
- 输入：所有基因按 score 排序（无需 P 值 cutoff）
- 方法：Kolmogorov-Smirnov 游走检验，评估基因集在排序列表顶部或底部是否富集
- 特点：捕获微弱但协同变化的基因集，无需人工设阈值
- 输出文件：`enrichment_gsea.csv` 及 `tables/enrichment/prerank_*_summary.csv`

### 主要输出文件

#### `tables/enrichment_ora.csv`（兼容入口）

汇总所有基因集的 ORA 结果。每个基因集对应 `tables/enrichment/ora_{gene_set}_summary.csv`。

| 列名 | 含义 |
|------|------|
| `cluster` | 细胞类型 |
| `Term` | 富集到的通路/GO 术语名称 |
| `Overlap` | 输入基因中命中该通路的比例（`命中数/通路总基因数`） |
| `Adjusted P-value` | 校正后 P 值（Fisher 精确检验 + Benjamini-Hochberg） |
| `Genes` | 命中基因列表 |
| `n_genes_input` | 该组输入的基因数量 |

#### `tables/enrichment_gsea.csv`（兼容入口）

汇总所有基因集的 pre-ranked GSEA 结果。每个基因集对应 `tables/enrichment/prerank_{gene_set}_summary.csv`。

| 列名 | 含义 |
|------|------|
| `cluster` | 细胞类型 |
| `Term` | 通路名称 |
| `NES` | 归一化富集分数（正=基因集在排序列表顶部富集，负=底部） |
| `NOM p-val` | 名义 P 值 |
| `FDR q-val` | 错误发现率校正后的 q 值 |
| `ES` | 原始富集分数 |
| `Leading edge` | 对富集贡献最大的核心基因 |

### 可视化图片

- `figures/enrichment/ora_{gene_set}_bubble.pdf` — ORA 气泡图
  - X 轴：细胞类型；Y 轴：通路名称
  - 点大小：Overlap（命中基因数）
  - 点颜色：`-log10(Adjusted P-value)`，越红越显著
- `figures/enrichment/prerank_{gene_set}_bubble.pdf` — GSEA 气泡图
  - X 轴：细胞类型；Y 轴：通路名称
  - 点大小：`-log10(FDR q-val)`
  - 点颜色：NES（红色=上调富集，蓝色=下调富集）

### 技术说明

- **网络要求**：首次运行时通过 Enrichr Web API 自动下载基因集库，缓存到本地（需要互联网）
- **基因集库**：支持 200+ 基因集库，常用包括：
  - `GO_Biological_Process` — 基因本体生物学过程
  - `KEGG_2021_Human` — KEGG 通路
  - `Reactome_2022` — Reactome 通路
  - `MSigDB_Hallmark_2020` — 标志性基因集
  - `WikiPathway_2021_Human` — WikiPathway
- **物种**：通过 `CFG.enrichment_organism` 配置（默认 Human），影响基因名匹配
- **缓存机制**：首次查询某基因集库后自动缓存，后续分析无需重复下载
- **容错设计**：当 `marker_genes_per_group.csv` 缺失时自动从 h5ad 重新计算标记基因再跑富集
- **并行 API 调用**：多个细胞类型的富集查询通过 `ThreadPoolExecutor` 并行发送，默认最多 5 线程
- **AI 解读**：如果启用 AI 配置，富集完成后会自动生成生物学解读报告保存为 `tables/enrichment/ai_interpretation.txt`

---

## 探索性分析输出

Step 11（`06_exploratory.py`）是一个可选的探索性步骤，生成以下辅助材料：

### 细胞组成分析

- `figures/06_exploratory/composition_by_stage_*.png` — 各类细胞随发育阶段的比例变化堆叠柱状图
- `tables/composition_by_stage_*.csv` — 对应的数值表

### QC 与标记基因可视化

- `figures/06_exploratory/_06_qc_umap.pdf` — UMAP 上显示 n_genes、total_counts、pct_counts_mt 等 QC 指标
- `figures/06_exploratory/_06_markers_*.pdf` — 已知标记基因在 UMAP 上的表达（每批最多 12 个基因）
- `figures/06_exploratory/_06_marker_dotplot.pdf` — 已知标记基因的 dotplot 概览

### 聚类统计

- `tables/cell_type_sizes.csv` — 各细胞类型的细胞数及占比
- `tables/leiden_sizes.csv` — 各 Leiden 聚类的细胞数及占比

---

## 输出文件速查表

| 分析类型 | 输出文件 | 生成步骤 |
|----------|----------|----------|
| 原始数据检查点 | `h5ad/00_raw.h5ad` | Step 00 |
| 双细胞检查点 | `h5ad/01_doublet.h5ad` | Step 02 |
| QC 检查点 | `h5ad/02_qc.h5ad` | Step 03 |
| 集成检查点 | `h5ad/03_integrated.h5ad` | Step 04 |
| 聚类检查点 | `h5ad/04_clustered.h5ad` | Step 05 |
| 注释检查点 | `h5ad/05_annotated.h5ad` | Step 06 |
| 轨迹检查点 | `h5ad/05_final.h5ad` | Step 09 |
| 标记基因 DE | `tables/marker_genes_per_group.csv` | Step 08 |
| 阶段配对 DE | `tables/pairwise_stage_de.csv` | Step 08 |
| 时间趋势 DE | `tables/temporal_trend_genes.csv` | Step 08 |
| 分支 DEG | `tables/branch_deg.csv` | Step 09 |
| ORA 富集 | `tables/enrichment_ora.csv` | Step 10 |
| GSEA 富集 | `tables/enrichment_gsea.csv` | Step 10 |
| 富集气泡图 | `figures/enrichment/` | Step 10 |
| 标记基因热图 | `figures/07_marker_heatmap.pdf` | Step 08 |
| 伪时间 UMAP | `figures/08_pseudotime.pdf` | Step 09 |
| 细胞组成 | `figures/06_exploratory/` | Step 11 |
