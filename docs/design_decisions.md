# 关键设计决策

> 关键设计决策详细参考 — 从 README 迁移至此

## 设计原则

- **集中配置 + CLI 主控 + h5ad 检查点链 + 断点恢复**：`config.py` 统一管理参数，`run_pipeline.py` 作为 CLI 入口控制执行流程，每一步写入独立 h5ad 文件形成检查点链，支持从任意中断步骤恢复。
- **HVG 在 raw counts 上选择**：`seurat_v3` 方法需要原始 counts，归一化会扭曲方差估计，必须在 raw 层计算高变基因。
- **保留 `.raw` 全基因表达**：即使 subset 到 HVG 进行降维，`.raw` 层仍保留所有基因的归一化表达，确保下游任意基因均可做差异表达和可视化。
- **每步只做一件事**：每个步骤脚本职责单一，输出的 h5ad 可独立加载、验证和调试。
- **不留技术债**：保留校正前后的 PCA 嵌入、不覆盖原始降维结果，便于对比和排查。

## 关键决策

### 分析方法

| 决策 | 说明 |
|------|------|
| HVG 在 raw counts 上选择 | 归一化会扭曲方差估计，`seurat_v3` 需要原始 counts |
| 保留 `.raw` | 确保下游任意基因可做 DE/可视化，不因 HVG subset 丢失信息 |
| Scrublet 双细胞检测 | per sample 并行，当前 scRNA-seq 最佳实践 |
| 复杂度指标 | `log10(genes) / log10(UMI)` 排除空滴/破损细胞 |
| 多参数网格聚类 | 6 分辨率 × 3 k 值，交互式选择最佳参数 |
| 分支轨迹分析 | PAGA 后比较相邻谱系的差异基因 |

### 工程架构

| 决策 | 说明 |
|------|------|
| 集中配置 + CLI | `config.py` + `run_pipeline.py` 避免硬编码，支持断点恢复 |
| 双模式 + 分块注释 | AI LLM 注释 + Score_genes 回退；大聚类数时分块提交防止超长 prompt |
| 富集回退 | 当 CSV 缺失时自动从 h5ad 计算标记基因再跑富集 |
| AI 交互工作流 | 5-Phase 交互式分析，agent 驱动 + 用户决策 |

### 实战淬炼

| 决策 | 说明 |
|------|------|
| 降采样 (Step 01) | 可选层化/随机/分样本封顶降采样，防止大数据集 OOM |
| NaN/Inf 数据完整性校验 | 各步骤间自动检测并修复 NaN/Inf，阻断传播 |
| `regress_out` 在 HVG 子集 + 顺序调优 | 先 normalize 再 `regress_out` 避免 NaN，子集上运行降内存 ~7× |
| 线粒体基因物种自定义 | `mt_gene_pattern` + `mt_gene_list` 双模式，支持非人/鼠物种 |
| Legacy 10X 兼容 | 自动检测 2 列 `genes.tsv.gz` 并补全 `feature_type` 列 |
| AI caller 重试 + reasoning 模型 | None 响应自动重试，支持 reasoning token 预算 |
| 并行计算优化 | ThreadPool/ProcessPool/joblib 并行化 UMAP、DE、富集 API 调用 |
