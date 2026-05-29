# 通用 scRNA-seq 分析管线

基于 GSE169109（人胎下丘脑）和 GSE138002（人视网膜）两个真实 scRNA-seq 项目的代码审查，提取出的一套通用、可复用的单细胞RNA测序分析工作流。

## 设计原则

- **工程架构学 GSE138002**：集中配置（config.py）、CLI 主控（run_pipeline.py）、h5ad checkpoint 链、支持断点恢复
- **分析方法学 GSE169109**：Scrublet 双细胞检测、raw counts 上选 HVG、子聚类、分支轨迹分析
- **每步只做一件事**，每一步的 h5ad 都可独立加载验证
- **不留技术债**：保留 `.raw` 全基因表达、不覆盖原始 PCA、同时保留校正前后嵌入

## 管线步骤

| 步骤 | 脚本 | 输入 | 输出 | 功能 |
|------|------|------|------|------|
| 00 | `00_load.py` | 原始数据 (MTX/CSV/h5ad) | `00_raw.h5ad` | 加载三种格式 + 样本/阶段映射 |
| 01 | `01_qc.py` | `00_raw.h5ad` | `01_qc.h5ad` | QC 指标 + Scrublet + 过滤 |
| 02 | `02_normalize_hvg.py` | `01_qc.h5ad` | `02_normalized.h5ad` | 归一化 + HVG + 保存 `.raw` |
| 03 | `03_pca_harmony.py` | `02_normalized.h5ad` | `03_harmony.h5ad` | PCA + Harmony 批次校正 |
| 04 | `04_cluster_umap.py` | `03_harmony.h5ad` | `04_clustered.h5ad` | 邻居图 + UMAP + 多分辨率 Leiden |
| 05 | `05_annotate.py` | `04_clustered.h5ad` | (同文件更新) | 标记基因自动注释 + 子聚类 |
| 06 | `06_exploratory.py` | `04_clustered.h5ad` | CSV + 图片 | 组成分析 + QC + 标记探索 |
| 07 | `07_markers_de.py` | `04_clustered.h5ad` | CSV + 图片 | 标记基因 + 配对DE + 时间趋势 |
| 08 | `08_trajectory.py` | `04_clustered.h5ad` | `05_final.h5ad` | PAGA + DPT + 分支分析 |
| 09 | `09_enrichment.py` | `marker_genes_per_group.csv` | `enrichment_*.csv` + 图片 | GO/KEGG 通路富集 (GSEApy) |

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
CFG.data_format = '10X_mtx'          # '10X_mtx' | 'csv_matrix' | 'h5ad'
CFG.mtx_prefix = 'GSE169109_'        # MTX 文件前缀
CFG.sample_map = {1: 'sample1', ...} # barcode 后缀 → 样本名
CFG.stage_map = {1: 'GW7', ...}      # barcode 后缀 → 发育阶段
CFG.stage_order = ['GW7', 'GW8', ...]
CFG.marker_dict = {                   # 已知标记基因
    'NE': ['HES1', 'SOX2', 'NES'],
    'Neuron': ['RBFOX3', 'MAP2', 'DCX'],
    ...
}
```

### 3. 运行

```bash
# 全部执行
python run_pipeline.py --config config_myproject.py

# 从断点恢复
python run_pipeline.py --resume --config config_myproject.py

# 只跑第 3 步
python run_pipeline.py --step 3 --config config_myproject.py

# 跑 3-6 步
python run_pipeline.py --steps 3-6 --config config_myproject.py
```

### 4. 输出结构

```
results/
├── h5ad/
│   ├── 00_raw.h5ad
│   ├── 01_qc.h5ad
│   ├── 02_normalized.h5ad
│   ├── 03_harmony.h5ad
│   ├── 04_clustered.h5ad
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

## 关键设计决策（相对于原项目）

| 决策 | 说明 | 来源项目 |
|------|------|----------|
| HVG 在 raw counts 上选择 | 归一化会扭曲方差估计，seurat_v3 需要原始 counts | GSE169109 |
| 保留 `.raw` | 确保下游任意基因可做 DE/可视化，不因 HVG subset 丢失信息 | GSE169109 |
| Scrublet 双细胞检测 | per sample 并行，是当前 scRNA-seq 最佳实践 | GSE169109 |
| 复杂度指标 | log10(genes)/log10(UMI) 排除空滴/破损细胞 | GSE169109 |
| 多分辨率 Leiden | 6 个分辨率 (0.3~2.0)，比 3 个更健壮 | GSE138002 |
| 集中配置 + CLI | config.py + run_pipeline.py 避免硬编码，支持断点恢复 | GSE138002 |
| 分支轨迹分析 | PAGA 后比较相邻谱系的差异基因 | GSE169109 |
| 时间趋势 Spearman | 连续发育阶段的相关性分析，不限于离散比较 | GSE169109 |
| 富集分析 | GSEApy (Enrichr API) + pre-ranked GSEA，双模式覆盖 | 新增 (GSEApy) |
| 富集回退 | 当 CSV 缺失时自动从 h5ad 计算标记基因再跑富集 | 新增 |

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

## 适配新数据集

1. 复制 `config.py` → `config_newdata.py`
2. 设置数据路径和格式
3. 配置 `sample_map`（如果是 10X 聚合数据）或 `meta_columns`（如果是 CSV 格式）
4. 配置 `marker_dict`（目标组织的已知细胞类型标记基因）
5. 调整 QC 阈值（建议先用默认值跑一遍看 QC 指标分布再微调）
6. 运行 `python run_pipeline.py --config config_newdata.py`

## 参考

- GSE169109：人胎下丘脑发育 scRNA-seq（Kim et al.）
- GSE138002：人视网膜发育 scRNA-seq（Sridhar et al.）
- Scanpy：https://scanpy.readthedocs.io/
- Harmony (PyTorch)：https://github.com/lilab-bcb/harmony-pytorch
- Scrublet：https://github.com/AllonKleinLab/scrublet
- GSEApy：https://gseapy.readthedocs.io/ (Bioinformatics 2022)
- Enrichr：https://maayanlab.cloud/Enrichr/
