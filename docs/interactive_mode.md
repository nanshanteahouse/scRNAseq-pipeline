# 交互模式参考文档

> 交互式工作流详细参考 — 从 README 迁移至此

## 概述

交互模式通过 `SKILL.md`（注册为 OpenCode skill `scRNAseq-interactive`）提供，是一套 AI Agent 驱动的 5-Phase 交互式 scRNA-seq 分析工作流。AI Agent 充当分析助手，在每个关键决策节点停下来与用户讨论，由用户确认后再继续。计算步骤自动通过 `run_pipeline.py` 执行，生物学解读由 AI 辅助完成。

## 工作流全景

```
Phase 0: 初始化（数据配置 + 资源规划：CPU 核数/OOM 恢复/降采样策略）
    │
    ▼
Phase 1: 数据预处理（自动，交互式降采样）
  downsample (可选) → doublet → QC → integrate
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

## 核心交互模式

```
AI 跑计算 → 展示结果 + 提出问题 → 用户决策 → AI 执行决策 → 继续下一段
```

## Phase 详情

| Phase | 名称 | 描述 | 用户决策点 |
|-------|------|------|------------|
| 0 | 初始化 | 数据路径/格式确认、组织类型、实验设计、资源规划（CPU核数、降采样策略、OOM恢复） | 数据配置确认、资源参数选择 |
| 1 | 数据预处理 | 降采样（可选）→ Scrublet → QC → Integrate | STOP 1: QC阈值决策 |
| 2 | 降维+聚类 | 多参数网格搜索（n_neighbors × resolution） | STOP 2: 最佳参数选择 |
| 3 | Major注释 | AI双模式注释 + Score_genes回退 | STOP 3: 注释结果审核/修改 |
| 4 | 亚型分析 | 对选定细胞类型循环执行：重聚类+AI重注释 | STOP 4: 亚型确认+DEG设计 |
| 5 | 下游分析 | DEG → 轨迹 → 富集 → AI生物学解读 | STOP 5: 结果报告审核 |

## 如何运行

```bash
python run_pipeline.py --interactive --config config_myproject.py
```

在 OpenCode 中加载 skill 后，agent 会自动按 5-Phase 交互工作流驱动分析，无需手动逐步骤操作。
