#!/usr/bin/env python3
"""
config.py — 通用 scRNA-seq 分析管线集中配置
=============================================

设计原则:
  - 所有参数集中在一个 dataclass 中，避免硬编码散落在各脚本里
  - 继承自 GSE169109（下丘脑）和 GSE138002（视网膜）两项目的最佳实践
  - 支持 10X MTX / CSV 矩阵 / h5ad 三种输入格式
  - 12 步管线 (Step 00-10): raw → QC → doublet → norm → pca/harmony → integrate → cluster → annotate → DE → trajectory → enrichment → report
  - 路径自动解析: 默认所有路径相对于 config.py 所在目录

使用方法:
    from config import CFG
    CFG.resolve_paths()
    adata = sc.read(CFG.norm_h5ad)

自定义配置:
    复制本文件为 config_myproject.py，修改参数后:
    python run_pipeline.py --config config_myproject.py
"""

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class AIConfig:
    """AI configuration — all AI features controlled here, disabled by default"""
    # Global switch
    enabled: bool = False

    # Inference endpoint (api_base is the sole switching mechanism)
    #   "http://<local_lan_url>"       → local vLLM (default, free)
    #   "https://api.deepseek.com/v1"  → DeepSeek API
    #   "http://localhost:11434/v1"    → Ollama
    api_base: str = ""
    model: str = "deepseek-chat"
    api_key: str = ""
    max_tokens: int = 4096
    temperature: float = 0.1

    # Task-level switches
    ai_qc_review: bool = False
    ai_param_suggest: bool = False
    ai_annotation: bool = True
    ai_subcluster: bool = True
    ai_deg_design: bool = False
    ai_interpretation: bool = True


@dataclass
class Config:
    # ═══════════════════════════════════════════════════════════════════
    #  路径设置
    # ═══════════════════════════════════════════════════════════════════
    data_dir: str = ""                   # 数据目录（自动设为 config.py 所在目录）
    results_dir: str = "results"         # 输出根目录
    h5ad_dir: str = "results/h5ad"       # AnnData checkpoint 目录
    figure_dir: str = "results/figures"  # 图片输出目录
    table_dir: str = "results/tables"    # CSV 表格输出目录
    log_dir: str = "logs"               # 日志文件目录

    # ═══════════════════════════════════════════════════════════════════
    #  数据输入格式
    # ═══════════════════════════════════════════════════════════════════
    # 可选: '10X_mtx'  |  'csv_matrix'  |  'h5ad'
    data_format: str = "10X_mtx"

    # 10X MTX 格式 (data_format='10X_mtx')
    mtx_prefix: str = ""                 # 文件前缀，如 "GSE169109_"
    mtx_dir: str = ""                    # MTX 文件目录（默认 DATA_DIR）

    # CSV 矩阵格式 (data_format='csv_matrix')
    matrix_file: str = ""                # 计数矩阵 (.mtx 或 .csv)
    barcodes_file: str = ""              # 细胞条形码文件
    features_file: str = ""              # 基因名文件

    # h5ad 格式 (data_format='h5ad')
    input_h5ad: str = ""                 # 直接读取已有的 h5ad
    backed: str = ""                     # h5ad 延迟加载模式: ''=全量, 'r'=只读backed (超大数据集用)

    # ═══════════════════════════════════════════════════════════════════
    #  样本元数据映射
    # ═══════════════════════════════════════════════════════════════════
    # barcode 后缀 → 样本名 → 发育阶段的映射
    # 10X多通道聚合时，barcode 以 -1, -2, ... 后缀区分样本
    # GSE169109 示例: {1: 'GW7-lane1', 2: 'GW7-lane2', ...}
    sample_map: Dict[int, str] = field(default_factory=dict)

    # barcode 后缀 → 发育阶段
    stage_map: Dict[int, str] = field(default_factory=dict)

    # 发育阶段的有序列表（用于图例排序和时间趋势分析）
    stage_order: List[str] = field(default_factory=list)

    # CSV 导入时细胞元数据中已有的列名映射（仅 data_format='csv_matrix'）
    # 例如: {'sample_col': 'sample', 'stage_col': 'age', 'tissue_col': 'sample_type'}
    meta_columns: Dict[str, str] = field(default_factory=dict)

    # ═══════════════════════════════════════════════════════════════════
    #  QC 阈值
    # ═══════════════════════════════════════════════════════════════════
    min_genes: int = 500                 # 细胞最少检测到的基因数
    max_genes: int = 7500                # 细胞最多检测到的基因数（排除可能双细胞）
    max_pct_mito: float = 20.0           # 线粒体基因百分比上限
    min_genes_per_umi: float = 0.7       # 复杂度下限 log10(genes)/log10(UMI)
    min_cells_per_gene: int = 3          # 基因至少在多少个细胞中表达

    # 自适应 MAD 阈值（GSE169109 策略）
    # 开启后，按每个样本的 median ± n*MAD 计算阈值，更适应异质性组织
    use_adaptive_thresholds: bool = False
    mad_n_mads: float = 3.0

    # ═══════════════════════════════════════════════════════════════════
    #  Scrublet 双细胞检测
    # ═══════════════════════════════════════════════════════════════════
    run_scrublet: bool = True
    scrublet_expected_doublet_rate: float = 0.06
    scrublet_min_counts: int = 2
    scrublet_min_cells: int = 3
    scrublet_min_gene_var_pctl: int = 85
    scrublet_n_prin_comps: int = 30

    # ═══════════════════════════════════════════════════════════════════
    #  归一化与 HVG
    # ═══════════════════════════════════════════════════════════════════
    normalize_target_sum: float = 1e4    # 总 counts 归一化目标值
    n_top_genes: int = 4000              # HVG 数量
    hvg_flavor: str = "seurat_v3"        # 'seurat_v3' | 'seurat' | 'cell_ranger'
    hvg_batch_key: str = "sample"        # 批次感知 HVG 的分组列名

    # ═══════════════════════════════════════════════════════════════════
    #  PCA
    # ═══════════════════════════════════════════════════════════════════
    n_pcs_full: int = 100                # 先计算 100 PC 用于 elbow 判断
    n_pcs_use: int = 50                  # 下游实际使用的 PC 数

    # ═══════════════════════════════════════════════════════════════════
    #  Harmony 批次校正
    # ═══════════════════════════════════════════════════════════════════
    use_harmony: bool = True
    harmony_batch_key: str = "sample"    # 批次协变量列名
    harmony_max_iter: int = 20

    # ═══════════════════════════════════════════════════════════════════
    #  聚类与 UMAP
    # ═══════════════════════════════════════════════════════════════════
    n_neighbors: int = 30                # 邻居图 k 值
    leiden_resolutions: List[float] = field(
        default_factory=lambda: [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
    )
    # Parameter grid for multi-run clustering (Step 04)
    param_grid_n_neighbors: list = field(default_factory=lambda: [15, 20, 30])
    param_grid_resolutions: list = field(default_factory=lambda: [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
    leiden_flavor: str = "igraph"        # 'igraph' | 'leidenalg'
    best_resolution: float = 1.0         # 注释时使用的主分辨率

    # ═══════════════════════════════════════════════════════════════════
    #  细胞类型注释
    # ═══════════════════════════════════════════════════════════════════
    # 格式: {'CellType': ['marker1', 'marker2', ...]}
    marker_dict: Dict[str, List[str]] = field(default_factory=dict)

    # 需要子聚类的细胞类型（如 ['NE', 'NP', 'Neuron']）
    subcluster_types: List[str] = field(default_factory=list)
    subcluster_resolution: float = 0.4
    min_cells_subcluster: int = 50

    # ═══════════════════════════════════════════════════════════════════
    #  差异表达分析
    # ═══════════════════════════════════════════════════════════════════
    de_method: str = "wilcoxon"          # 'wilcoxon' | 't-test'
    de_n_genes: int = 50                 # 每组 top N 基因
    de_pval_cutoff: float = 0.05
    de_logfc_cutoff: float = 0.25        # 设为 None 则不过滤 logFC

    # ═══════════════════════════════════════════════════════════════════
    #  轨迹分析
    # ═══════════════════════════════════════════════════════════════════
    # 根细胞类型列表（如 ['NE', 'NP']）
    root_cell_types: List[str] = field(default_factory=list)
    # 自动检测根时使用的标记基因（如 ['SOX2', 'PAX6', 'HES1']）
    root_markers: List[str] = field(default_factory=list)
    n_diffmap_comps: int = 15
    n_branchings: int = 2

    # ═══════════════════════════════════════════════════════════════════
    #  富集分析 (Step 09)
    # ═══════════════════════════════════════════════════════════════════
    run_enrichment: bool = True            # 是否运行富集分析
    # 富集方法: 'ora' (over-representation) | 'prerank' (pre-ranked GSEA) | 'both'
    enrichment_method: str = "both"
    # 基因集库: 见 https://maayanlab.cloud/Enrichr/#libraries
    # 常用: 'GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
    #       'KEGG_2021_Human', 'Reactome_2022', 'MSigDB_Hallmark_2020',
    #       'WikiPathway_2023_Human'
    enrichment_gene_sets: list = field(
        default_factory=lambda: [
            'GO_Biological_Process_2023',
            'KEGG_2021_Human',
        ]
    )
    enrichment_organism: str = "human"     # 'human' | 'mouse' | 'rat' ...
    enrichment_n_top_genes: int = 200      # ORA 模式: 每个聚类取 top N 上调基因
    enrichment_pval_cutoff: float = 0.05   # 富集结果 p-value 阈值
    enrichment_min_size: int = 10          # 基因集最小大小
    enrichment_max_size: int = 500         # 基因集最大大小
    enrichment_permutations: int = 1000    # GSEA 置换次数

    # ═══════════════════════════════════════════════════════════════════
    #  执行环境
    # ═══════════════════════════════════════════════════════════════════
    n_jobs: int = 24
    random_seed: int = 42
    scanpy_verbosity: int = 2            # 0=quiet, 1=warn, 2=info, 3=hint
    force_csr: bool = True               # 强制 X 为 CSR 格式（行优先，细胞级操作最优）
    use_float32: bool = False            # 强制 float32 精度（节省 ~50% 内存，默认关闭）

    # ═══════════════════════════════════════════════════════════════════
    #  AI 配置
    # ═══════════════════════════════════════════════════════════════════
    ai: AIConfig = field(default_factory=AIConfig)

    # ═══════════════════════════════════════════════════════════════════
    #  中间 h5ad checkpoint 路径（自动派生）
    # ═══════════════════════════════════════════════════════════════════
    @property
    def raw_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "00_raw.h5ad")

    @property
    def qc_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "02_qc.h5ad")

    @property
    def doublet_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "01_doublet.h5ad")

    @property
    def norm_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "02_normalized.h5ad")

    @property
    def harmony_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "03_harmony.h5ad")

    @property
    def integrated_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "03_integrated.h5ad")

    @property
    def cluster_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "04_clustered.h5ad")

    @property
    def annotated_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "05_annotated.h5ad")

    @property
    def final_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "05_final.h5ad")

    # ──────────────────────────────────────────────────────────────────
    #  方法
    # ──────────────────────────────────────────────────────────────────
    def resolve_paths(self):
        """
        解析所有路径。
        非绝对路径将被视为相对于调用者（config.py）所在目录。
        """
        base = os.path.dirname(os.path.abspath(__file__))
        for attr in [
            "data_dir", "results_dir", "h5ad_dir",
            "figure_dir", "table_dir", "log_dir",
            "mtx_dir",
        ]:
            val = getattr(self, attr)
            if val and not os.path.isabs(val):
                setattr(self, attr, os.path.join(base, val))

        # 默认 data_dir = config.py 所在目录
        if not self.data_dir:
            self.data_dir = base
        if not self.mtx_dir:
            self.mtx_dir = self.data_dir

        # 创建输出目录
        for d in [self.results_dir, self.h5ad_dir,
                  self.figure_dir, self.table_dir, self.log_dir]:
            os.makedirs(d, exist_ok=True)

    def has_sample_mapping(self) -> bool:
        """是否配置了 barcode suffix → sample 映射"""
        return len(self.sample_map) > 0

    def has_stage_mapping(self) -> bool:
        """是否配置了 stage 映射"""
        return len(self.stage_map) > 0

    def has_markers(self) -> bool:
        """是否配置了标记基因"""
        return len(self.marker_dict) > 0


# ═══════════════════════════════════════════════════════════════════════
#  全局实例
# ═══════════════════════════════════════════════════════════════════════
# 各脚本统一通过 from config import CFG 导入
CFG = Config()
