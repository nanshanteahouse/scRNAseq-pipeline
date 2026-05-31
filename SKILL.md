# SKILL: scRNA-seq 交互式分析管线

> 供 AI Agent 加载，用于在对话中驱动 scRNA-seq 全流程分析。
> 核心原则：**计算是流水线，解读才是价值。工作流 = 计算 run ↔ 交互 stop 交替进行。**

---

## 技能概述

本技能定义了一套 **5-Phase 交互式 scRNA-seq 分析工作流**。AI Agent 充当分析助手，在每个关键决策节点停下来与用户讨论，由用户确认后再继续。计算步骤自动执行，生物学解读由 AI 辅助完成。

### 你（AI Agent）的角色

- **分析驱动者** —— 依次推进 5 个 Phase，每次跑完一个计算块就停在交互节点
- **生物学顾问** —— 用 LLM（本地 vLLM / DeepSeek API）辅助注释、解读结果
- **对话引导者** —— 在交互节点主动展示结果、提出问题、给出建议，等待用户决策

### 用户交互模式

```
AI 跑计算 → 展示结果 + 提出问题 → 用户决策 → AI 执行决策 → 继续下一段
```

---

## 工作流全景

```
00_load
    │
    ▼
┌────────────────────────────────────────────────────┐
│  Phase 1: 数据预处理（自动跑，无交互）                │
│  01_doublet → 02_qc → 03_integrate                  │
└────────────────────────────────────────────────────┘
    │
    ▼
────────────────  INTERACTIVE STOP 1 ────────────────
  QC 审查：展示 doublet 率、过滤统计 → 用户决策
────────────────────────────────────────────────────
    │
    ▼
┌────────────────────────────────────────────────────┐
│  Phase 2: 降维 + 聚类（多参数网格）                   │
│  04_pca → 05_multi_cluster                          │
└────────────────────────────────────────────────────┘
    │
    ▼
────────────────  INTERACTIVE STOP 2 ────────────────
  Layout 选择：展示各参数 UMAP → AI 推荐 → 用户选择
────────────────────────────────────────────────────
    │
    ▼
┌────────────────────────────────────────────────────┐
│  Phase 3: Major type 注释（AI 协作）                  │
│  06_annotate_major                                   │
└────────────────────────────────────────────────────┘
    │
    ▼
────────────────  INTERACTIVE STOP 3 ────────────────
  注释审核：展示 AI 标注结果 → 用户确认/修改
────────────────────────────────────────────────────
    │
    ▼
┌────────────────────────────────────────────────────┐
│  Phase 4: 亚型分析（循环，AI 协作）                    │
│  FOR EACH selected type: 07_subcluster              │
└────────────────────────────────────────────────────┘
    │
    ▼
────────────────  INTERACTIVE STOP 4 ────────────────
  亚型确认 + DEG 方案讨论
────────────────────────────────────────────────────
    │
    ▼
┌────────────────────────────────────────────────────┐
│  Phase 5: 下游分析 + AI 解读                          │
│  08_deg → 09_trajectory → 10_enrichment              │
└────────────────────────────────────────────────────┘
    │
    ▼
────────────────  INTERACTIVE STOP 5 ────────────────
  结果解读：展示 DEG/GO/轨迹 → AI 生成生物学解读报告
────────────────────────────────────────────────────
```

---

## 步骤注册表

```python
STEPS = [
    ("00", "00_load.py",           "加载原始数据 → 00_raw.h5ad"),
    ("01", "01_doublet.py",        "Scrublet 双细胞检测 (per sample) → 01_doublet.h5ad"),
    ("02", "02_qc.py",             "QC 过滤 (先去 doublet) → 02_qc.h5ad"),
    ("03", "03_integrate.py",      "归一化 + HVG + PCA + Harmony → 03_integrated.h5ad"),
    ("04", "04_pca.py",            "PCA 降维 → 04_pca.h5ad"),
    ("05", "05_multi_cluster.py",  "多参数 UMAP + 多分辨率 Leiden → checkpoint 待确认"),
    ("06", "06_annotate_major.py", "AI 辅助 major cell type 注释"),
    ("07", "07_subcluster.py",     "交互式亚型分析 (循环)"),
    ("08", "08_deg.py",            "差异表达分析 (多层级)"),
    ("09", "09_trajectory.py",     "PAGA + DPT 轨迹分析"),
    ("10", "10_enrichment.py",     "GO/KEGG 富集 + AI 解读"),
]
```

---

## Phase 0：初始化

在开始任何分析前，确认以下信息：

**必须获取：**
1. **数据路径** — 原始数据在哪里？（matrix/MTX/h5ad）
2. **数据格式** — `10X_mtx` / `csv_matrix` / `h5ad`？
3. **组织 & 物种** — 用于 AI 注释（如 "human brain", "mouse retina"）
4. **实验设计** — 有哪些分组/条件？用户关心的比较是什么？

**可选获取：**
5. **样本映射** — barcode 后缀 → 样本名
6. **阶段/时间点** — 如果是发育数据
7. **已知 marker 基因** — 如果有先验知识

**初始化操作：**
```python
# 复制 config.py 为 config_project.py → 修改配置
# 检查 00_load.py 能否正确读取数据
python run_pipeline.py --step 00 --config config_project.py
```

---

## Phase 1：数据预处理（自动）

此阶段不需要交互，一次性跑完。

### Step 01 — doublet（Per sample Scrublet）

```bash
python run_pipeline.py --step 01 --config config_project.py
```

- 对每个样本独立运行 Scrublet
- 输出 `01_doublet.h5ad`，含 `doublet_scores`, `predicted_doublet` 列

### Step 02 — QC 过滤

```bash
python run_pipeline.py --step 02 --config config_project.py
```

- 先去除标记为 doublet 的细胞
- 再按 `min_genes` / `max_genes` / `max_pct_mito` / `min_genes_per_umi` 过滤
- 输出 `02_qc.h5ad`

### Step 03 — 归一化 + HVG + PCA + Harmony

```bash
python run_pipeline.py --step 03 --config config_project.py
```

- `normalize_total` → `log1p` → 选 HVG → `regress_out` → `scale` → `pca` → `harmony`
- 输出 `03_integrated.h5ad`

---

## ⏸️ INTERACTIVE STOP 1 — QC 审查

### 你做的事情

1. 自动读取 01 和 02 步的输出日志/结果文件
2. 向用户展示以下信息：

```
📊 QC 完成，概况如下：

Sample        Total cells  Doublets(%)  After filter
sample_A      12,340       823 (6.7%)    10,112
sample_B      15,678       1,021 (6.5%)  12,890
sample_C      8,901        356 (4.0%)    7,234
─────────────────────────────────────
Total         36,919       2,200 (6.0%)  30,236

MT% 分布: 中位 8.2%，99th pctl 22.1%
基因数分布: 中位 2,300，<500 的细胞占 3.1%
复杂度分布: 中位 0.82，<0.7 的细胞占 5.4%
```

3. **询问用户**：

```
QC 阈值需要调整吗？当前：
  • min_genes=500  max_genes=7500
  • max_pct_mito=20%
  • min_genes_per_umi=0.7

你可以：
  ✅ "继续" — 接受当前阈值
  ✏️ "把 MT 阈值改成 15" — 调整后 rerun
  ❓ "MT% 偏高正常吗？" — 我解释
```

### 等待用户决策

- 用户说"继续" → 锁定 checkpoint，进入 Phase 2
- 用户调整阈值 → 修改 config → rerun Step 02 → 回到此 stop
- 用户追问 → 用 LLM 生成解释（调用 `ai_qc_review` prompt）

---

## Phase 2：降维 + 聚类（半自动）

### Step 04 — PCA

```bash
python run_pipeline.py --step 04 --config config_project.py
```

- 计算 100 PC
- 输出 elbow plot
- 输出 `04_pca.h5ad`

### Step 05 — 多参数网格聚类

```bash
python run_pipeline.py --step 05 --config config_project.py
```

- 跑一组参数网格：

```python
param_grid = {
    "n_neighbors": [15, 20, 30],
    "leiden_resolution": [0.3, 0.5, 0.8, 1.0, 1.5, 2.0],
}
```

- 每组参数输出 UMAP 图 + 聚类统计
- **不写 checkpoint**（等用户选了再写）

---

## ⏸️ INTERACTIVE STOP 2 — Layout 选择

### 你做的事情

1. 收集所有参数组合的 UMAP 图片 + 聚类统计
2. 调用 LLM（`ai_param_suggest` prompt）让 AI 推荐最佳参数
3. 向用户展示：

```
🗺️ 聚类结果 — 建议你重点看这几组：

┌──────────────┬──────────┬──────────┬──────────────────┐
│ n_neighbors  │ resolution│ n_clusters │ silhouette_score │
├──────────────┼──────────┼──────────┼──────────────────┤
│ 20           │ 0.8      │ 12       │ 0.52 ⭐最高       │
│ 15           │ 1.0      │ 16       │ 0.48             │
│ 30           │ 0.5      │ 8        │ 0.45             │
│ 20           │ 1.5      │ 22       │ 0.40             │
└──────────────┴──────────┴──────────┴──────────────────┘

🤖 AI 推荐：n_neighbors=20, resolution=0.8
  → 12 个 cluster，silhouette 最高，分离度好
  → [UMAP 图]

参数网格 UMAP 对比：
  [图1: res=0.5]  [图2: res=0.8]  [图3: res=1.0]  [图4: res=1.5]
```

4. **询问用户**：

```
请选择你要用的参数：
  • "用 0.8" — 锁定 resolution=0.8
  • "n=30, res=1.0 看看" — 我要那组结果
  • "我想试一个新的参数" — 你说我跑
  • "继续" — 用 AI 推荐的默认
```

### 等待用户决策

- 用户选择 → 锁定 checkpoint，保存选定的参数到 config
- 用户想试新参数 → 补充跑一步 → 回到此 stop
- **只有确认后才写 `04_clustered.h5ad`**

---

## Phase 3：Major Type 注释（AI 协作）

### Step 06 — AI 注释

```bash
# 自动触发 AI 注释，或：
python run_pipeline.py --step 06 --config config_project.py
```

**AI 注释流程：**

```
1. 读 04_clustered.h5ad
2. sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
3. 取每群 top 20 marker 基因
4. 构造 prompt → 发送给 LLM
5. 解析返回的 JSON → 写入 adata.obs
6. 输出 marker_genes.csv + 候选注释
```

**AI 不直接写入**，结果先缓存供用户审核。

---

## ⏸️ INTERACTIVE STOP 3 — 注释审核

### 你做的事情

1. 展示 AI 注释结果 + 推理依据

```
🏷️ AI 注释结果：

Cluster   Cell Type        State         Confidence   Key Markers
───────   ─────────        ─────         ──────────   ──────────
0         T cell           activated     🔴 high     CD3D, CD2, CD69
1         Microglia        resting       🔴 high     CSF1R, CD74, P2RY12
2         Oligodendrocyte  N/A           🟡 medium   MBP, PLP1, SOX10
3         T cell           N/A           🟢 low      (mix of T & Myeloid markers)
4         Astrocyte        N/A           🔴 high     GFAP, AQP4, SLC1A2
5         Neuron           excitatory    🟡 medium   SLC17A7, SATB2
...

⚠️ 低置信度 cluster:
  • Cluster 3 (low) — marker 不明确，~30% 细胞混合
```

2. **询问用户**：

```
你可以：
  ✅ "确认" — 用 AI 的结果
  ✏️ "Cluster 3 改成 'T cell activated'" — 手动修改
  ❓ "为什么 Cluster 3 分不好？" — 我解释
  🔄 "重新跑 Cluster 3 的注释" — 调整 prompt 再试
```

### 等待用户决策

- 用户确认 → 锁定注释 → 写回 `processed_adata.h5ad`
- 用户修改个别 cluster → 应用修改 → 展示更新
- 用户想重跑 → 调整 prompt（如"专注 T cell marker"）→ 重新调用 AI

---

## Phase 4：亚型分析（循环）

### Step 07 — 交互式子聚类

对用户选定的每个 major type，执行以下循环：

```python
for cell_type in selected_types:
    # 1. 提取子集
    sub = adata[adata.obs['cell_type'] == cell_type]
    
    # 2. 重新 PCA + neighbors + multi-resolution UMAP + Leiden
    sc.pp.neighbors(sub, n_pcs=50)
    sc.tl.leiden(sub, resolution=[0.2, 0.4, 0.8])
    
    # 3. AI 重新注释亚型
    markers = sc.get.rank_genes_groups_df(sub, group=None)
    # 调用 LLM → 返回 subtype 注释
    
    # 4. 停在交互节点让用户确认
    # → 用户确认后写回 adata.obs['cell_type_sub']
```

---

## ⏸️ INTERACTIVE STOP 4 — 亚型确认 + DEG 方案

### 对每个 cell type：

```
🧬 子聚类: T cell (3,240 cells)

可选分辨率:
  res=0.2 → 3 亚型  [UMAP]
  res=0.4 → 5 亚型  [UMAP] ⭐
  res=0.8 → 8 亚型  [UMAP]

🤖 AI 推荐 res=0.4:
  • Subtype 0: CD8+ T_cytotoxic (CD8A, GZMB, PRF1)
  • Subtype 1: CD4+ T_naive (CD4, CCR7, LEF1, SELL)
  • Subtype 2: T_reg (FOXP3, IL2RA, CTLA4)
  • Subtype 3: CD8+ T_exhausted (PDCD1, LAG3, TOX)
  • Subtype 4: NK (NKG7, KLRD1, GNLY)

确认？(yes/modify/skip)
```

**所有 major type 完成后，讨论 DEG 方案：**

```
📋 DEG 实验设计

adata.obs 中可用分组列:
  • condition: [WT, KO]
  • treatment: [Control, Drug]
  • sample: [S1, S2, S3, S4]

推荐方案:
  • 两两比较: KO vs WT（每个 cell type 内）
  • 或线性模型: ~ condition * treatment

你要怎么比较？
```

---

## Phase 5：下游分析

### Step 08 — DEG

```bash
python run_pipeline.py --step 08 --config config_project.py
```

- 对 major type + subtype 各跑一次
- 输出 `deg_results_cell_type.csv`, `deg_results_cell_subtype.csv`
- 火山图 per group comparison

### Step 09 — 轨迹分析

```bash
python run_pipeline.py --step 09 --config config_project.py
```

- PAGA + DPT + 分支 DEG + 表达趋势

### Step 10 — 富集分析 + AI 解读

```bash
python run_pipeline.py --step 10 --config config_project.py
```

- Enrichr API (GO BP/MF/CC, KEGG, Reactome)
- 热图 + 柱状图
- 最终：**调用 LLM 生成生物学解读报告**

---

## ⏸️ INTERACTIVE STOP 5 — 结果解读

### 你做的事情

1. 收集 DEG / GO / 轨迹的关键结果
2. 调用 LLM（`ai_interpretation` prompt）生成跨细胞类型的解读
3. 向用户展示解读报告

```
🧬 生物学解读报告 — 项目: hypoxia_mouse_brain

### 关键发现
1. **Microglia** 在缺氧后显著激活
   - 上调: IL1B, NLRP3, CCL2, CXCL10 (p.adj<0.001)
   - 富集: "Inflammatory Response", "NLRP3 Inflammasome"
   - 解读: 缺氧诱导的神经炎症主要通过 Microglia 介导

2. **Oligodendrocyte** 显示髓鞘修复相关变化
   - 上调: PLP1, MBP, CLDN11 (log2FC≈1.2-1.8)
   - 富集: "Myelination", "Axon Ensathment"
   - 解读: 缺氧后少突胶质细胞启动髓鞘再生程序

3. **T cell** 浸润增加
   - 在缺氧组中 T cell 比例从 2% 增至 8%
   - 上调: CD8A, GZMB, IFNG — 细胞毒性 T 细胞活化

### 可验证假设
- Microglia 中的 NLRP3 通路是缺氧后神经炎症的核心驱动
- 阻断 CCL2-CCR2 轴可能减少 T 细胞浸润
```

4. **询问用户**：

```
你可以：
  ✅ "完成" — 保存报告
  ❓ "T cell 那个有文献支持吗？" — 我补充
  🔄 "只看 Microglia 的结果，换 GO 数据库再跑一次" — rerun
```

---

## AI 集成参考

### 配置

```python
# config.py 中
@dataclass
class AIConfig:
    enabled: bool = False
    api_base: str = "http://<local_lan_url>"  # 或 api.deepseek.com/v1
    model: str = "DeepSeek-V4-Flash"
    api_key: str = ""
    max_tokens: int = 4096
    temperature: float = 0.1
    # 任务级开关
    ai_annotation: bool = True
    ai_subcluster: bool = True
    ai_interpretation: bool = True
```

### 调用

```python
from scripts.ai_caller import ai_query

response = ai_query(
    system_prompt="You are an expert scRNA-seq biologist...",
    user_prompt="Cluster markers: ...\nTissue: brain\nSpecies: human",
    cfg=cfg.ai,
)
```

### Prompt 模板

见 `scripts/ai_prompts.py`，包含：
- `ANNOTATION_SYSTEM_PROMPT` / `ANNOTATION_USER_PROMPT_TEMPLATE` — 细胞注释
- `PARAM_SUGGEST_PROMPT` — 参数推荐（待实现）
- `QC_REVIEW_PROMPT` — QC 审查（待实现）
- `DEG_DESIGN_PROMPT` — 实验设计（待实现）
- `INTERPRETATION_PROMPT` — 生物学解读（待实现）

---

## 错误处理

### 计算步骤失败

1. 读取错误日志（`logs/*.log`）
2. 向用户展示错误概要
3. 如果 AI 修复开启（`ai_auto_heal=True`），调用 LLM 诊断并尝试修复
4. 否则询问用户是否要重试或跳过

```
❌ Step 03 (integrate) 失败

错误: OSError — 内存不足，无法分配 8GB array

建议:
  • 增加 `n_jobs=8` 减少并行
  • 或者用 `n_pcs_full=50` 减少 PC 数
  • 或者在有更多内存的节点上运行

要调整后重试吗？
```

### LLM 调用失败

1. 检查 `api_base` 是否可达
2. 如果本地服务挂了，询问用户是否 fallback 到 DeepSeek API
3. 如果所有 LLM 不可用，回退到 `score_genes` 方法

---

## 快速开始

```bash
# 1. 配数据
cp config.py config_myproject.py
# 编辑 data_dir, data_format, sample_map, ...

# 2. 跑加载
python run_pipeline.py --step 00 --config config_myproject.py

# 3. 以交互模式启动完整管线
python run_pipeline.py --interactive --config config_myproject.py
```

交互模式下，AI Agent（即你）负责在每次 stop 节点展示结果、引导讨论、等待决策。
