# 新数据集适配指南

> 适配新数据集详细参考 — 从 README 迁移至此

---

## 快速检查清单

1. **复制配置文件** — `cp config.py config_newdata.py`
2. **设置数据格式** — 确定原始数据格式，修改 `CFG.data_format`（`10X_mtx` / `csv_matrix` / `h5ad`）
3. **配置样本映射** — 10X 聚合数据配置 `CFG.sample_map`（barcode 后缀 → 样本名）；CSV 格式配置 `CFG.meta_columns`
4. **配置发育阶段** — 设置 `CFG.stage_map`（barcode 后缀 → 发育阶段）和 `CFG.stage_order`（有序列表）
5. **配置标记基因** — 在 `CFG.marker_dict` 中填写目标组织的已知细胞类型标记基因
6. **调整 QC 阈值** — 建议先用默认值跑一遍，观察 QC 指标分布后再微调
7. **考虑降采样** — 大数据集（>50K 细胞）设置 `CFG.downsample_target` 激活降采样；小数据集跳过（默认）
8. **运行管线** — `python run_pipeline.py --config config_newdata.py`

---

## 数据格式详解

管线支持三种输入格式，通过 `CFG.data_format` 切换。

### 10X_mtx（标准 Cell Ranger 输出）

```
mtx_dir/
├── myproject_matrix.mtx.gz        # 计数矩阵
├── myproject_features.tsv.gz      # 基因特征（3 列：id, symbol, feature_type）
└── myproject_barcodes.tsv.gz      # 细胞条形码
```

- `CFG.mtx_prefix`：文件前缀，如 `myproject_`
- `CFG.mtx_dir`：MTX 文件所在目录（默认等于 `CFG.data_dir`）
- 加载后自动从 barcode 后缀（如 `-1`, `-2`）解析样本和阶段映射
- 支持 **legacy 旧版格式**：自动检测 2 列 `genes.tsv.gz` 并补全 `feature_type` 列（加 `\tGene Expression`），无需手动处理
- 清理 `gene_ids` 列，避免下游干扰

### csv_matrix（自定义 CSV 格式）

适用于非 10X 标准格式的数据，如公共数据集中常见的 CSV 矩阵 + 元数据文件。

**两种子模式：**

| 子模式 | 触发条件 | 说明 |
|--------|----------|------|
| True CSV | `matrix_file` 后缀为 `.csv` 或 `.gz` | 基因 × 细胞的 CSV 表，第一列为基因名。加载后自动转置为细胞 × 基因 |
| MTX + 元数据 | `matrix_file` 后缀为 `.mtx` 或 `.mtx.gz` | 用 `mmread` 加载 MTX 矩阵，通过 `barcodes_file` 读入元数据 |

核心参数：

- `CFG.matrix_file`：计数矩阵文件路径
- `CFG.barcodes_file`：细胞元数据 CSV（`index_col=0`），可包含 sample / stage / tissue 等列
- `CFG.features_file`：基因名列表（第一列作为基因名）
- `CFG.meta_columns`：列名映射字典，将 CSV 中的列重命名为管线标准名（`sample`, `stage`, `tissue`）

**欧洲逗号处理**：部分公共数据集的元数据使用欧洲格式逗号（`,`）作为小数分隔符（如 `12,5` 表示 12.5）。管线通过 pandas CSV 解析器自动兼容处理。

### h5ad（预处理的 AnnData）

直接加载已有 h5ad 文件，适合：
- 之前处理过的数据，想复用下游分析步骤
- 从其他管线导出为 h5ad 格式的数据

核心参数：

- `CFG.input_h5ad`：h5ad 文件路径
- `CFG.backed`：可选，设为 `'r'` 使用 AnnData backed 模式，适合超大数据集（不全部加载到内存）

> 注意：h5ad 格式不会自动添加样本/阶段映射，如需映射需在加载脚本中手动处理或提前在 h5ad 中包含这些列。

---

## 配置字段速查

以下列出适配新数据集时必须关注的关键配置字段：

```python
# ── 数据格式 ──
CFG.data_format = '10X_mtx'          # 原始数据格式：10X_mtx | csv_matrix | h5ad
CFG.mtx_prefix = 'myproject_'         # MTX 文件前缀（10X_mtx 格式）

# ── 样本与阶段映射 ──
CFG.sample_map = {1: 'sample1', ...} # barcode 后缀 → 样本名（10X 聚合数据）
CFG.stage_map = {1: 'GW7', ...}      # barcode 后缀 → 发育阶段
CFG.stage_order = ['GW7', 'GW8', ...] # 阶段排序（图例和趋势分析用）
CFG.meta_columns = {                  # CSV 格式的元数据列映射
    'sample': 'sample',
    'stage': 'age',
    'tissue': 'sample_type',
}

# ── 先验知识 ──
CFG.tissue = 'hypothalamus'          # 组织类型
CFG.species = 'human'                # 物种（影响线粒体基因检测）
CFG.marker_dict = {                   # 已知细胞类型标记基因
    'NE': ['HES1', 'SOX2', 'NES'],
    'Neuron': ['RBFOX3', 'MAP2', 'DCX'],
    ...
}

# ── QC 阈值（默认值适合大部分数据集） ──
CFG.min_genes = 500                  # 细胞最少基因数
CFG.max_genes = 7500                 # 细胞最多基因数
CFG.max_pct_mito = 20.0             # 线粒体基因百分比上限
CFG.min_genes_per_umi = 0.7         # 复杂度下限

# ── 聚类 ──
CFG.leiden_resolutions = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
CFG.best_resolution = 1.0           # 注释使用的主分辨率

# ── 轨迹分析 ──
CFG.root_cell_types = ['NE', 'NP']  # 根细胞类型
CFG.root_markers = ['SOX2', 'PAX6'] # 自动检测根时使用的标记基因

# ── 执行环境 ──
CFG.n_jobs = 24                     # 并行计算核数
```

---

## QC 阈值调优建议

默认阈值在大多数数据集中表现良好，但最优值因组织、物种和质量而异。

**推荐做法：**

1. 先用默认阈值跑一次 `run_pipeline.py --steps 0-3` 完成加载到 QC 过滤
2. 检查 `02_qc.h5ad` 的 QC 指标分布：
   ```python
   import scanpy as sc
   adata = sc.read('results/h5ad/02_qc.h5ad')
   sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])
   ```
3. 根据实际分布调整阈值：
   - **`min_genes`**：观察 n_genes_by_counts 的左尾，确保不切掉真实细胞
   - **`max_genes`**：观察右尾，排除可能的双细胞
   - **`max_pct_mito`**：线粒体比例过高通常表示细胞破损，但某些组织（如肝脏）本底高
   - **`min_genes_per_umi`**：复杂度低可能表示空滴或技术噪声
4. 更新阈值后重新跑步骤 2-3（管线支持断点恢复）

> 也可开启 `CFG.use_adaptive_thresholds = True`，使用每样本 median ± n*MAD 的自适应阈值策略，更适合异质性组织。

---

## 降采样策略

### 何时启用

- **数据集 >50K 细胞**：建议降采样以加速实验和节省内存
- **数据集 >200K 细胞**：强烈建议降采样，否则后续步骤（UMAP、Leiden、轨迹）可能 OOM

### 如何配置

```python
CFG.downsample_target = 50000       # 目标总细胞数
CFG.downsample_strategy = 'stratified'  # 层化降采样
```

设置 `CFG.downsample_target` 后 Step 01（`downsample.py`）自动激活；不设置则 Step 01 自动跳过（默认行为）。

### 三种策略

| 策略 | 取值 | 说明 | 适用场景 |
|------|------|------|----------|
| 随机 | `random` | 从全体细胞中均匀随机采样 | 各类样本比例均衡时 |
| 层化 | `stratified` | 按样本分层采样，保持各样本比例 | 样本间细胞数差异大时（推荐） |
| 分样本封顶 | `max_per_sample` | 每样本最多取 N 个细胞，其余丢弃 | 大样本压迫小样本时 |

```python
# max_per_sample 需要额外参数
CFG.downsample_strategy = 'max_per_sample'
CFG.downsample_max_per_sample = 10000  # 每样本上限
```

### 注意事项

- 降采样在 Step 01 执行，**直接覆写 `00_raw.h5ad`**，下游步骤无感知
- 开启后如需完整数据，请先备份 `00_raw.h5ad`
- 降采样随机种子固定（`CFG.downsample_random_seed = 42`），保证可复现

---

## 特殊场景处理

### 非人物种

默认线粒体基因检测使用 `MT-` 前缀模式（适用于人类）。对非人/鼠物种需要调整：

```python
# 模式匹配法：使用正则前缀
CFG.mt_gene_pattern = 'mt-'          # 小鼠线粒体基因前缀（小写）

# 或列表法：明确指定线粒体基因（优先于 pattern）
CFG.mt_gene_list = ['ATP6', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
```

`mt_gene_list` 优先级高于 `mt_gene_pattern`，同时设置时以列表为准。

### 旧版 10X 数据（2 列 genes.tsv.gz）

10X 早期版本输出 2 列 `genes.tsv.gz`（gene_id + symbol），缺少 `feature_type` 列。管线 00_load.py 自动检测此类文件：

```python
# 自动执行：检测到 genes.tsv.gz 且无 features.tsv.gz 时
# 复制 genes.tsv.gz 为 features.tsv.gz，每行末尾追加 '\tGene Expression'
```

用户无需任何额外配置，管线自动处理。

### 欧洲逗号分隔符

部分公共数据的元数据 CSV 使用逗号 `,` 作为小数分隔符（如 `12,5` 而非 `12.5`）。通用管线通过 pandas CSV 解析器自动兼容处理。

---

## 交互模式说明

在交互式模式下（`run_pipeline.py --interactive`），资源决策由 AI Agent 在 Phase 0 自动引导用户完成：

- CPU 核数配置
- 降采样策略选择
- OOM 恢复方案

上述决策无需在配置文件中手动设置，Agent 会根据数据集大小和用户环境自动推荐。

---

## 模板配置文件参考

项目中包含两个实际项目的配置模板，可直接参考：

| 模板文件 | 数据格式 | 关键特点 |
|----------|----------|----------|
| `templates/config_templates/config_hypothalamus.py` | 10X_mtx | 多通道聚合数据，多个发育阶段 |
| `templates/config_templates/config_retina.py` | csv_matrix | CSV 格式 + 元数据列映射 |

使用方式：
```bash
cp templates/config_templates/config_hypothalamus.py config_myproject.py
# 然后修改其中的配置参数
```
