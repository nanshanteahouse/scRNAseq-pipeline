# 通用 scRNA-seq 分析管线

基于多个真实 scRNA-seq 项目的工程与方法论提炼，覆盖**加载 → QC → 降采样 → 聚类 → 注释 → DE → 轨迹 → 富集**全流程的开放管线。

**两种运行模式：**

- **Pipeline 模式** — 全自动串行执行，CLI 控制，支持断点恢复，适合批处理
- **Interactive 模式** — AI Agent 驱动的 5-Phase 交互工作流，适合探索性分析

---

## 设计原则

- **集中配置 + CLI 主控**：`config.py` + `run_pipeline.py`，h5ad checkpoint 链，支持断点恢复
- **每步只做一件事**：每步的 h5ad 可独立加载验证
- **不留技术债**：保留 `.raw` 全基因表达、校正前后嵌入、不覆盖原始 PCA
- **AI Agent 原生支持**：预定义 5-Phase 交互工作流 + LLM 细胞注释 + score_genes 双模式回退

---

## 快速开始

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

> **WSL 用户注意**：数据在 `/mnt/` 路径下需设置 `export HDF5_USE_FILE_LOCKING=FALSE` 避免 h5py 文件锁定。

### 2. 配置

复制 `config.py` 为 `config_myproject.py`，设置以下关键参数：

```python
# 数据输入
CFG.data_format = '10X_mtx'          # '10X_mtx' | 'csv_matrix' | 'h5ad'
CFG.mtx_prefix = 'myproject_'        # MTX 文件前缀

# 样本元数据
CFG.sample_map = {1: 'sample1', ...} # barcode 后缀 → 样本名
CFG.stage_map = {1: 'GW7', ...}      # barcode 后缀 → 发育阶段
CFG.stage_order = ['GW7', 'GW8', ...]

# 先验知识
CFG.marker_dict = {                   # 已知细胞类型标记基因
    'NE': ['HES1', 'SOX2', 'NES'],
    'Neuron': ['RBFOX3', 'MAP2', 'DCX'],
}
```

参考模板：`templates/config_templates/config_hypothalamus.py`、`templates/config_templates/config_retina.py`

### 3. 运行

**Pipeline 模式**：

```bash
# 全流程执行
python run_pipeline.py --config config_myproject.py

# 从断点恢复
python run_pipeline.py --resume --config config_myproject.py

# 单步执行 (step 4 = 04_integrate.py)
python run_pipeline.py --step 4 --config config_myproject.py

# 范围执行
python run_pipeline.py --steps 4-7 --config config_myproject.py

# 列出所有步骤
python run_pipeline.py --list --config config_myproject.py
```

**Interactive 模式**：

```bash
python run_pipeline.py --interactive --config config_myproject.py
```

详见 [交互式工作流文档](docs/interactive_mode.md)。

### 4. 输出结构

```
results/
├── h5ad/         # checkpoint 链
├── figures/      # 图表输出
└── tables/       # CSV 结果表
```

各文件详细说明见 [输出文件参考](docs/output_reference.md)。

---

## 项目文件结构

```
.
├── run_pipeline.py          # CLI 主控
├── config.py                # 集中配置
├── SKILL.md                 # AI Agent skill（交互工作流）
├── opencode.json            # OpenCode skill 注册
├── requirements.txt         # 依赖
├── scripts/                 # 管线脚本
│   ├── 00_load.py
│   ├── downsample.py        # step 01：降采样（可选）
│   ├── 01_doublet.py
│   ├── 02_qc.py
│   ├── 03_integrate.py
│   ├── 04_cluster_umap.py
│   ├── 05_annotate_major.py # AI + Score_genes 双模式
│   ├── 06_subcluster.py
│   ├── 06_exploratory.py    # step 11（可选）
│   ├── 07_markers_de.py
│   ├── 08_trajectory.py
│   ├── 09_enrichment.py
│   ├── ai_caller.py         # LLM 调用封装
│   ├── ai_prompts.py        # Prompt 模板
│   └── utils.py             # 通用工具函数
├── projects/                # 项目文件夹
├── templates/
│   └── config_templates/    # 配置模板
├── docs/                    # 详细文档
│   ├── pipeline_steps.md
│   ├── interactive_mode.md
│   ├── design_decisions.md
│   ├── ai_integration.md
│   ├── new_dataset.md
│   └── output_reference.md
├── logs/                    # 日志（已 gitignore）
├── results/                 # 输出（已 gitignore）
└── cache/                   # 缓存（已 gitignore）
```

---

## 详细文档

| 文档 | 内容 |
|------|------|
| [管线步骤详解](docs/pipeline_steps.md) | 12 步脚本说明、checkpoint 链、CLI 命令、关键注意事项 |
| [交互式工作流](docs/interactive_mode.md) | 5-Phase AI Agent 驱动分析流程、各阶段决策点 |
| [关键设计决策](docs/design_decisions.md) | 分析方法/工程架构/实战淬炼各层面的设计选择与理由 |
| [AI 集成配置](docs/ai_integration.md) | AIConfig 详解、后端支持、双模式注释、回退机制 |
| [适配新数据集](docs/new_dataset.md) | 数据格式、配置步骤、QC 调优、特殊场景处理 |
| [输出文件说明](docs/output_reference.md) | DE 分层、富集结果、checkpoint 文件结构 |

---

## 参考

- [Scanpy](https://scanpy.readthedocs.io/)
- [Harmony (PyTorch)](https://github.com/lilab-bcb/harmony-pytorch)
- [Scrublet](https://github.com/AllonKleinLab/scrublet)
- [GSEApy](https://gseapy.readthedocs.io/)
- [Enrichr](https://maayanlab.cloud/Enrichr/)
