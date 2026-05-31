# 通用 scRNA-seq 分析管线

基于 GSE169109（人胎下丘脑）和 GSE138002（人视网膜）两个真实 scRNA-seq 项目的代码审查，提取出的一套通用、可复用的单细胞 RNA 测序分析工作流。

**两种模式：**
- **Pipeline 模式** — 全自动串行执行，适合批处理场景
- **Interactive 模式** — AI Agent 驱动的 5-Phase 交互工作流，适合探索性分析

---

## 设计原则

- **工程架构学 GSE138002**：集中配置（`config.py`）、CLI 主控（`run_pipeline.py`）、h5ad checkpoint 链、支持断点恢复
- **分析方法学 GSE169109**：Scrublet 双细胞检测、raw counts 上选 HVG、子聚类、分支轨迹分析
- **每步只做一件事**，每一步的 h5ad 都可独立加载验证
- **不留技术债**：保留 `.raw` 全基因表达、不覆盖原始 PCA、同时保留校正前后嵌入
- **AI Agent 原生支持**：预定义 5-Phase 交互式工作流，agent 自动驱动分析 + 用户协作决策

---

## 管线步骤

| 步骤 | 脚本 | 输入 | 输出 | 功能 |
|------|------|------|------|------|
| 00 | `00_load.py` | 原始数据 (MTX/CSV/h5ad) | `00_raw.h5ad` | 加载三种格式 + 样本/阶段映射 |
| 01 | `01_doublet.py` | `00_raw.h5ad` | `01_doublet.h5ad` | Scrublet 双细胞检测 (per sample) |
| 02 | `01_qc.py` | `01_doublet.h5ad` | `01_qc.h5ad` | QC 过滤 (先去 doublet，再按 MT%/基因数/复杂度过滤) |
| 03 | `03_integrate.py` | `01_qc.h5ad` | `03_integrated.h5ad` | 归一化 + HVG + PCA + Harmony 批次校正 |
| 04 | `04_cluster_umap.py` | `03_integrated.h5ad` | `04_clustered.h5ad` | 多参数 UMAP + 多分辨率 Leiden 聚类 |
| 05 | `05_annotate_major.py` | `04_clustered.h5ad` | `05_annotated.h5ad` | AI/Score_genes 双模式 major type 注释 |
| 06 | `06_subcluster.py` | `05_annotated.h5ad` | (同文件更新) | 交互式亚型分析 (需指定 cell type) |
| 07 | `07_markers_de.py` | `05_annotated.h5ad` | CSV + 图片 | 标记基因 + 配对DE + 时间趋势 |
| 08 | `08_trajectory.py` | `04_clustered.h5ad` | `05_final.h5ad` | PAGA + DPT + 分支分析 |
| 09 | `09_enrichment.py` | `marker_genes_per_group.csv` | `enrichment_*.csv` + 图片 | GO/KEGG 通路富集 (GSEApy) |
| 10 | `06_exploratory.py` | `05_annotated.h5ad` | CSV + 图片 | 组成分析 + QC + 标记探索 |

> **关于 step 02**：步骤索引为 `2`，但脚本复用了 `01_qc.py`（doublet 检测 + QC 过滤作为独立的两步共用同一入口，通过配置控制行为）。

---

## Interactive 模式（AI Agent 驱动）

通过 `SKILL.md`（注册为 OpenCode skill `scRNAseq-interactive`）提供 AI Agent 驱动的 5-Phase 交互式分析工作流。

工作流全景：

```
Phase 0: 初始化（确认数据路径/格式/组织/实验设计）
    │
    ▼
Phase 1: 数据预处理（自动）
  doublet → QC → integrate
    │
    ▼
  STOP 1 — QC 审查：展示过滤统计 → 用户决策阈值
    │
    ▼
Phase 2: 降维 + 聚类（半自动）
  多参数 grid: n_neighbors x resolution
    │
    ▼
  STOP 2 — Layout 选择：展示各参数 UMAP → AI 推荐 → 用户选择
    │
    ▼
Phase 3: Major type 注释（AI 协作）
  AI 双模式注释 → 结果缓存
    │
    ▼
  STOP 3 — 注释审核：展示 AI 标注 → 用户确认/修改
    │
    ▼
Phase 4: 亚型分析（AI 协作，循环）
  FOR EACH selected type: re-cluster + AI re-annotation
    │
    ▼
  STOP 4 — 亚型确认 + DEG 实验设计
    │
    ▼
Phase 5: 下游分析 + AI 解读
  DEG → Trajectory → Enrichment
    │
    ▼
  STOP 5 — 结果解读：AI 生成生物学解读报告
```

核心交互模式：`AI 跑计算 → 展示结果 + 提出问题 → 用户决策 → AI 执行决策 → 继续下一段`

---

## 快速开始

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

> **WSL 用户注意**：如果数据在 `/mnt/` 路径下，需要设置环境变量避免 h5py 文件锁定：
> ```bash
> export HDF5_USE_FILE_LOCKING=FALSE
> ```

### 2. 配置

复制 `config.py` 为 `config_myproject.py`，修改以下关键参数：

```python
# 数据输入
CFG.data_format = '10X_mtx'          # '10X_mtx' | 'csv_matrix' | 'h5ad'
CFG.mtx_prefix = 'GSE169109_'        # MTX 文件前缀

# 样本元数据
CFG.sample_map = {1: 'sample1', ...} # barcode 后缀 → 样本名
CFG.stage_map = {1: 'GW7', ...}      # barcode 后缀 → 发育阶段
CFG.stage_order = ['GW7', 'GW8', ...]

# 先验知识
CFG.marker_dict = {                   # 已知标记基因
    'NE': ['HES1', 'SOX2', 'NES'],
    'Neuron': ['RBFOX3', 'MAP2', 'DCX'],
    ...
}
```

> 项目中已包含模板配置供参考：`templates/config_templates/config_hypothalamus.py`、`templates/config_templates/config_retina.py`

### 3. 运行

**Pipeline 模式**（全自动执行）：

```bash
# 全部执行
python run_pipeline.py --config config_myproject.py

# 从断点恢复
python run_pipeline.py --resume --config config_myproject.py

# 只跑第 3 步
python run_pipeline.py --step 3 --config config_myproject.py

# 跑 3-6 步
python run_pipeline.py --steps 3-6 --config config_myproject.py

# 列出所有步骤
python run_pipeline.py --list --config config_myproject.py
```

**Interactive 模式**（AI Agent 驱动）：

在 OpenCode 中加载 skill 后，agent 会自动按 5-Phase 交互工作流驱动分析。详见 `SKILL.md`。

```bash
python run_pipeline.py --interactive --config config_myproject.py
```

### 4. 输出结构

```
results/
├── h5ad/
│   ├── 00_raw.h5ad
│   ├── 01_doublet.h5ad
│   ├── 01_qc.h5ad
│   ├── 03_integrated.h5ad
│   ├── 04_clustered.h5ad
│   ├── 05_annotated.h5ad
│   └── 05_final.h5ad
├── figures/
│   ├── pca_elbow.png
│   ├── harmony_comparison.png
│   ├── umap_leiden_resolutions.pdf
│   ├── 05_celltype.pdf
│   ├── 06_exploratory/
│   ├── 07_marker_heatmap.pdf
│   ├── 08_pseudotime.pdf
│   └── 09_enrichment/
│       ├── enrichment_*_ORA.png
│       └── enrichment_*_GSEA.png
└── tables/
    ├── marker_genes_per_group.csv
    ├── pairwise_stage_de.csv
    ├── temporal_trend_genes.csv
    ├── branch_deg.csv
    ├── enrichment_ora.csv          ← Step 09 ORA 结果
    └── enrichment_gsea.csv         ← Step 09 pre-ranked GSEA 结果
```

---

## 关键设计决策

| 决策 | 说明 | 来源 |
|------|------|------|
| HVG 在 raw counts 上选择 | 归一化会扭曲方差估计，seurat_v3 需要原始 counts | GSE169109 |
| 保留 `.raw` | 确保下游任意基因可做 DE/可视化，不因 HVG subset 丢失信息 | GSE169109 |
| Scrublet 双细胞检测 | per sample 并行，当前 scRNA-seq 最佳实践 | GSE169109 |
| 复杂度指标 | log10(genes)/log10(UMI) 排除空滴/破损细胞 | GSE169109 |
| 多参数网格聚类 | 6 分辨率 x 3 k 值，交互式选择最佳参数 | GSE138002 |
| 集中配置 + CLI | config.py + run_pipeline.py 避免硬编码，支持断点恢复 | GSE138002 |
| 分支轨迹分析 | PAGA 后比较相邻谱系的差异基因 | GSE169109 |
| 时间趋势 Spearman | 连续发育阶段的相关性分析，不限于离散比较 | GSE169109 |
| 双模式注释 | AI LLM 注释 + Score_genes 回退，确保离线也能用 | 新增 |
| 富集分析 | GSEApy (Enrichr API) + pre-ranked GSEA，双模式覆盖 | 新增 |
| 富集回退 | 当 CSV 缺失时自动从 h5ad 计算标记基因再跑富集 | 新增 |
| AI 交互工作流 | 5-Phase 交互式分析，agent 驱动 + 用户决策 | 新增 |

---

## 输出说明

- **Layer 1 DE** (`marker_genes_per_group.csv`)：每类细胞 vs 其他所有细胞的标记基因
- **Layer 2 DE** (`pairwise_stage_de.csv`)：同一类细胞在相邻发育阶段的差异基因
- **Layer 3 DE** (`temporal_trend_genes.csv`)：表达量随发育时间单调上升/下降的基因
- **Branch DEG** (`branch_deg.csv`)：轨迹分支决策点的谱系特化基因
- **ORA 富集** (`enrichment_ora.csv`)：每类上调标记基因的 GO/KEGG 过表达分析（超几何检验）
- **GSEA 富集** (`enrichment_gsea.csv`)：每类全基因排序的 GSEA（无需 cutoff，捕获协同变化）

> 富集分析通过 GSEApy 调用 Enrichr Web API（需网络）。首次分析会自动缓存基因集库。
> 支持 200+ 基因集库：GO_Biological_Process、KEGG、Reactome、MSigDB Hallmarks 等。
> 结果包含 `Adjusted P-value`（ORA）和 `FDR q-val` / `NES`（GSEA）供下游筛选。

---

## AI 集成

管线集成 AI 能力用于细胞注释、亚型分析和结果解读，通过 `config.py` 中的 `AIConfig` 控制。

```python
@dataclass
class AIConfig:
    enabled: bool = False
    api_base: str = ""               # 本地 vLLM / DeepSeek API / Ollama
    model: str = "deepseek-chat"
    api_key: str = ""
    # 任务级开关
    ai_annotation: bool = True       # 细胞注释
    ai_subcluster: bool = True       # 亚型分析
    ai_interpretation: bool = True   # 结果解读
```

支持多种后端：
- **本地 vLLM**：`api_base = "http://<local_ip>:8000/v1"`
- **DeepSeek API**：`api_base = "https://api.deepseek.com/v1"`
- **Ollama**：`api_base = "http://localhost:11434/v1"`

---

## 适配新数据集

1. 复制 `config.py` → `config_newdata.py`
2. 设置数据路径和格式
3. 配置 `sample_map`（如果是 10X 聚合数据）或 `meta_columns`（如果是 CSV 格式）
4. 配置 `marker_dict`（目标组织的已知细胞类型标记基因）
5. 调整 QC 阈值（建议先用默认值跑一遍看 QC 指标分布再微调）
6. 运行 `python run_pipeline.py --config config_newdata.py`

---

## 项目文件结构

```
.
├── run_pipeline.py          # CLI 主控
├── config.py                # 集中配置
├── SKILL.md                 # AI Agent skill（交互工作流）
├── opencode.json             # OpenCode skill 注册
├── requirements.txt         # 依赖
├── .gitignore
├── scripts/
│   ├── 00_load.py
│   ├── 01_doublet.py
│   ├── 01_qc.py             # 复用：step 01 doublet + step 02 QC
│   ├── 03_integrate.py
│   ├── 04_cluster_umap.py
│   ├── 05_annotate_major.py # AI + Score_genes 双模式
│   ├── 06_subcluster.py
│   ├── 06_exploratory.py    # step 10（可选）
│   ├── 07_markers_de.py
│   ├── 08_trajectory.py
│   ├── 09_enrichment.py
│   ├── ai_caller.py         # LLM 调用封装
│   ├── ai_prompts.py        # Prompt 模板
│   └── utils.py             # 通用工具函数
├── templates/
│   └── config_templates/    # 项目配置模板
│       ├── config_hypothalamus.py
│       └── config_retina.py
├── logs/                    # 日志（已 gitignore）
├── results/                 # 输出（已 gitignore）
├── cache/                   # 缓存（已 gitignore）
└── figures/                 # 图片输出
```

---

## 参考

- GSE169109：人胎下丘脑发育 scRNA-seq（Kim et al.）
- GSE138002：人视网膜发育 scRNA-seq（Sridhar et al.）
- Scanpy：https://scanpy.readthedocs.io/
- Harmony (PyTorch)：https://github.com/lilab-bcb/harmony-pytorch
- Scrublet：https://github.com/AllonKleinLab/scrublet
- GSEApy：https://gseapy.readthedocs.io/ (Bioinformatics 2022)
- Enrichr：https://maayanlab.cloud/Enrichr/
