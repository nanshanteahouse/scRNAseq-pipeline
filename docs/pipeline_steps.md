# 管线步骤详细参考 — 从 README 迁移至此

## 快速命令

```bash
# 完整运行
python run_pipeline.py --config config_myproject.py

# 单步运行（0-based 索引，与 STEPS 列表直接对应）
python run_pipeline.py --step 4 --config config_myproject.py

# 范围运行
python run_pipeline.py --steps 3-6 --config config_myproject.py

# 从第一个未完成的 checkpoint 恢复
python run_pipeline.py --resume --config config_myproject.py

# 列出所有步骤及其索引
python run_pipeline.py --list --config config_myproject.py
```

## 步骤总览

| 步骤 | 脚本 | 输入 | 输出 | 功能 |
|------|------|------|------|------|
| 00 | `00_load.py` | 原始数据 (MTX/CSV/h5ad) | `00_raw.h5ad` | 加载三种格式 + 样本/阶段映射 |
| 01 | `downsample.py` | `00_raw.h5ad` | (覆写 `00_raw.h5ad`) | 降采样 (可选，按配置跳过) |
| 02 | `01_doublet.py` | `00_raw.h5ad` | `01_doublet.h5ad` | Scrublet 双细胞检测 (per sample) |
| 03 | `02_qc.py` | `01_doublet.h5ad` | `02_qc.h5ad` | QC 过滤 (先去 doublet，再按 MT%/基因数/复杂度过滤) |
| 04 | `03_integrate.py` | `02_qc.h5ad` | `03_integrated.h5ad` | 归一化 + HVG + PCA + Harmony 批次校正 |
| 05 | `04_cluster_umap.py` | `03_integrated.h5ad` | `04_clustered.h5ad` | 多参数 UMAP + 多分辨率 Leiden 聚类 |
| 06 | `05_annotate_major.py` | `04_clustered.h5ad` | `05_annotated.h5ad` | AI/Score_genes 双模式 major type 注释 |
| 07 | `06_subcluster.py` | `05_annotated.h5ad` | (同文件更新) | 交互式亚型分析 (需指定 cell type) |
| 08 | `07_markers_de.py` | `05_annotated.h5ad` | CSV + 图片 | 标记基因 + 配对DE + 时间趋势 |
| 09 | `08_trajectory.py` | `04_clustered.h5ad` | `05_final.h5ad` | PAGA + DPT + 分支分析 |
| 10 | `09_enrichment.py` | `marker_genes_per_group.csv` | `enrichment_*.csv` + 图片 | GO/KEGG 通路富集 (GSEApy) |
| 11 | `06_exploratory.py` | `05_annotated.h5ad` | CSV + 图片 | 组成分析 + QC + 标记探索 |

> 所有脚本位于 `scripts/` 目录下。配置通过 `--config` 参数动态加载（`utils.resolve_config()`）。

## 步骤映射说明

CLI `--step N` 运行 `STEPS[N]`（0-based 索引）。例如 `--step 4` 运行 `04_integrate.py`。`--steps 3-6` 运行步骤 3、4、5、6。

| CLI 步骤 | 脚本 | 产出文件 | 读取文件 |
|----------|------|----------|----------|
| 0 | `00_load.py` | `00_raw.h5ad` | 原始数据文件 |
| 1 | `downsample.py` | (覆写 `00_raw.h5ad`) | `00_raw.h5ad` |
| 2 | `01_doublet.py` | `01_doublet.h5ad` | `00_raw.h5ad` |
| 3 | `02_qc.py` | `02_qc.h5ad` | `01_doublet.h5ad` |
| 4 | `03_integrate.py` | `03_integrated.h5ad` | `02_qc.h5ad` |
| 5 | `04_cluster_umap.py` | `04_clustered.h5ad` (通过 `04_grid_results.h5ad`) | `03_integrated.h5ad` |
| 6 | `05_annotate_major.py` | `05_annotated.h5ad` | `04_clustered.h5ad` |
| 7 | `06_subcluster.py` | (更新 `05_annotated.h5ad`) | `05_annotated.h5ad` |
| 8 | `07_markers_de.py` | `marker_genes_per_group.csv` | `05_annotated.h5ad` |
| 9 | `08_trajectory.py` | `05_final.h5ad` | `04_clustered.h5ad` |
| 10 | `09_enrichment.py` | `enrichment_*.csv` | `marker_genes_per_group.csv` |
| 11 | `06_exploratory.py` | CSV + 图片 | `05_annotated.h5ad` |

## 各步骤详解

### 步骤 00: 数据加载 (`00_load.py`)

加载三种格式的原始数据：
- **h5ad**: 直接读取。
- **MTX**: 10X 格式（`matrix.mtx`、`barcodes.tsv`、`features.tsv`），结合 CSV 元数据。
- **CSV**: 稠密矩阵。

自动应用 `sample_map` 和 `stage_map` 映射，为每个细胞标注样本来源和实验阶段。

### 步骤 01: 降采样 (`downsample.py`)

**可选步骤**。默认状态下自动跳过。

- 启用条件：在配置中设置 `CFG.downsample_target`（目标细胞总数）。
- 支持策略：`random`（随机）、`stratified`（分层）、`max_per_sample`（每样本最大数）。
- 启用后：读取 `00_raw.h5ad`，降采样后覆写同一文件。下游步骤（02+）透明地读取降采样后的数据。
- 禁用时：无操作，直接退出码 0。

### 步骤 02: 双细胞检测 (`01_doublet.py`)

对每个样本独立运行 Scrublet 进行双细胞（doublet）检测。使用 `joblib.Parallel` 并行处理多个样本。双细胞标签写入 `adata.obs['doublet']`。

### 步骤 03: QC 过滤 (`02_qc.py`)

按顺序执行以下过滤：
1. 去除标记为 doublet 的细胞。
2. 按线粒体基因百分比（MT%）过滤。
3. 按基因数（n_genes）过滤。
4. 按复杂度（复杂细胞 vs 空液滴/双细胞）过滤。

### 步骤 04: 批次校正与整合 (`03_integrate.py`)

完整流程：
1. 归一化（normalization）。
2. 高变基因（HVG）筛选。
3. **关键优化**：`regress_out` 在 HVG 子集上运行（约 4000 个基因），而非全基因组。此举将峰值内存降低约 7 倍。
4. PCA 降维。
5. Harmony 批次校正。

### 步骤 05: 聚类与 UMAP (`04_cluster_umap.py`)

多参数网格搜索：遍历多种 UMAP 参数和 Leiden 分辨率，将结果存入 `04_grid_results.h5ad`（中间文件）。

**注意**：管线预期的下游 checkpoint 是 `04_clustered.h5ad`，它在交互模式中由用户选择参数后生成。管线模式需固定分辨率。

### 步骤 06: Major Type 注释 (`05_annotate_major.py`)

双模式注释：
1. **AI 模式**：调用 OpenAI API 基于标记基因列表推断细胞类型。
2. **Score_genes 模式**（后备）：使用 `scanpy.tl.score_genes` 进行标记基因打分。

**AI 后备链**：优先尝试 AI 注释。遭遇连接失败、JSON 解析错误、键缺失等任何问题时，自动降级为 `score_genes` 模式。始终优雅处理，不会阻塞管线。

### 步骤 07: 亚型分析 (`06_subcluster.py`)

对指定细胞类型进行亚聚类分析。需要在配置中指定目标细胞类型。更新 `05_annotated.h5ad` 文件，新增亚型标签。

### 步骤 08: 标记基因与差异表达 (`07_markers_de.py`)

执行三种分析：
1. **标记基因检测**：识别各 cluster 的特异标记基因。
2. **配对差异表达**：组间对比分析。
3. **时间趋势分析**：沿时间/阶段检测表达变化。

输出：CSV 文件和可视化图片。

### 步骤 09: 轨迹分析 (`08_trajectory.py`)

基于 `04_clustered.h5ad` 运行：
- **PAGA**：基于 partition-based graph abstraction 的轨迹推断。
- **DPT**：扩散伪时间（diffusion pseudotime）。
- **分支分析**：识别分化分支点。

输出：`05_final.h5ad`。

### 步骤 10: 通路富集 (`09_enrichment.py`)

基于 `marker_genes_per_group.csv` 运行 GO/KEGG 通路富集分析。使用 GSEApy 库，通过 `ThreadPoolExecutor` 并行发起 API 请求。

**依赖注意**：gseapy > 0.11.0 需要 Rust 编译器。若 pip 安装失败，运行：

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

输出：`enrichment_*.csv` 和富集结果图片。

### 步骤 11: 探索性分析 (`06_exploratory.py`)

生成全面的探索性分析结果，包括：
- 细胞类型组成分析。
- QC 指标可视化。
- 标记基因表达探索。

输出：CSV 文件和多种图表。

## Checkpoint 链机制

所有 h5ad 文件构成一条线性的 checkpoint 链：

```
原始数据 → 00_raw.h5ad → 01_doublet.h5ad → 02_qc.h5ad → 03_integrated.h5ad → 04_clustered.h5ad → 05_annotated.h5ad
```

每个步骤读取上一个 checkpoint，处理后写入新的 checkpoint。恢复（resume）模式通过扫描步骤 0-6 的 checkpoint 文件，找到第一个缺失的从中断处继续。

### Safe Write 机制

`utils.safe_write()` 专为 WSL `/mnt` 文件锁定问题设计：
1. 先将 h5ad 写入 `/tmp` 目录下的临时文件。
2. 然后用 `mv` 命令原子性地移动到目标路径。

此举有效规避 WSL 下 h5py 文件锁定导致的写入失败。

## 关键注意事项

### 步骤编号偏移（v2 新增降采样步骤）

- CLI `--step N` 运行 `STEPS[N]`（0-based）。`--step 4` → `04_integrate.py`（原为 `--step 3`）。
- `--steps 4-7` 运行步骤 4、5、6、7。
- Checkpoint 文件名保持原编号：`01_doublet.h5ad` → `02_qc.h5ad` → `03_integrated.h5ad`。
- 配置属性：`CFG.qc_h5ad` = `02_qc.h5ad`，`CFG.doublet_h5ad` = `01_doublet.h5ad`。
- 以 `run_pipeline.py` 中的 `STEPS` 列表为权威映射。

### WSL h5py 文件锁定

如果数据位于 `/mnt/` 路径下，运行 Python **之前**必须设置：

```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

管线不会自动设置此环境变量。读取时若不设置可能失败；写入时 `safe_write()` 已做缓解，但读取仍需此设置。

### 步骤 04 输出混淆点

- `04_cluster_umap.py` 写入 `04_grid_results.h5ad`（中间结果，包含多参数网格搜索结果）。
- 下游步骤实际期望的 checkpoint 是 `04_clustered.h5ad`。
- 恢复（resume）逻辑检查 `05_annotated.h5ad` 来判断步骤 05 是否就绪，而非 `04_clustered.h5ad`。

### 步骤 06 AI 后备链

AI 注释在遇到以下情况时会自动降级：
- API 连接失败。
- JSON 解析错误。
- 返回数据中缺少必要字段。

降级后使用 `scanpy.tl.score_genes` 进行基于基因打分的传统注释，保证管线不中断。

### gseapy Rust 依赖

```bash
# 若 gseapy pip 安装失败，先安装 Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

## 日志系统

- `utils.setup_logger()` 同时输出到 stdout 和 `logs/<step>.log`。
- 日志格式：`%H:%M:%S | %-7s | %s`。

## 并行计算

| 步骤 | 并行方式 | 用途 |
|------|----------|------|
| 02 (Scrublet) | `joblib.Parallel` | 每样本独立检测 |
| 05 (聚类) | `ProcessPoolExecutor` | 多参数 UMAP/Leiden 网格搜索 |
| 10 (富集) | `ThreadPoolExecutor` | 并行 API 请求 |

## AI 集成

- `ai_caller.ai_query()` 封装 OpenAI SDK。
- 配置通过 `CFG.ai`（`AIConfig` 对象）。
- Prompt 定义在 `ai_prompts.py`。
- 目前仅步骤 05（注释）使用 AI 功能。
