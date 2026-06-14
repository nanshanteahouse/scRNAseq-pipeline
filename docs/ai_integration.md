# AI 集成参考

> AI 集成详细参考 — 从 README 迁移至此

## 概述

scRNA-seq 管线集成了 AI 能力，用于三大场景：

- **细胞类型注释**：基于 marker 基因表达自动推断聚类所属的细胞类型
- **亚型分析**：对指定细胞类型进行再聚类后重新注释亚型
- **结果解读**：对差异表达和富集分析结果生成自然语言生物学解读

所有 AI 功能通过 `config.py` 中的 `AIConfig` 控制，配置文件切换 `api_base` 即可切换后端。

---

## AIConfig 配置

```python
@dataclass
class AIConfig:
    enabled: bool = False
    api_base: str = ""               # 本地 vLLM / DeepSeek API / Ollama
    model: str = "deepseek-v4-flash"   # deepseek-chat 将于 2026/07/24 弃用
    api_key: str = ""
    # 任务级开关
    ai_annotation: bool = True       # 细胞注释
    ai_subcluster: bool = True       # 亚型分析
    ai_interpretation: bool = True   # 结果解读
```

各字段说明：

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `enabled` | `False` | 全局总开关。设为 `False` 时所有 AI 功能跳过，完全不调用 API |
| `api_base` | `""` | API 端点地址。**这也是切换后端的主要手段** — 改这个字段即可切换 vLLM / DeepSeek / Ollama |
| `model` | `"deepseek-v4-flash"` | 模型名称（deepseek-chat 将于 2026/07/24 弃用）。与 `api_base` 配合使用 |
| `api_key` | `""` | API 密钥。建议留空，通过 `.env` 文件或环境变量 `LLM_API_KEY` 提供（见 [安全存储](#api-key-安全存储)） |
| `max_tokens` | `4096` | 最大输出 token 数 |
| `temperature` | `0.1` | 采样温度（思考模式下不生效）。注释任务推荐低温（0.0–0.2） |
| `thinking_enabled` | `True` | 是否启用 DeepSeek 思考模式 |
| `reasoning_effort` | `"high"` | 推理强度：`"high"` 或 `"max"` |
| `timeout` | `None` | API 调用超时（秒）。`None` 表示不设超时 |
| `ai_annotation` | `True` | 是否启用 AI 细胞注释。关闭后回退到 `score_genes` 模式 |
| `ai_subcluster` | `True` | 是否启用 AI 亚型注释 |
| `ai_interpretation` | `True` | 是否启用 AI 结果解读 |
| `ai_qc_review` | `False` | 是否启用 AI QC 审查（预留字段） |
| `ai_param_suggest` | `False` | 是否启用 AI 参数建议（预留字段） |
| `ai_deg_design` | `False` | 是否启用 AI DEG 对比设计（预留字段） |

`AIConfig` 作为 `Config.ai` 属性挂载：

```python
CFG.ai.api_base = "http://localhost:11434/v1"
CFG.ai.model = "qwen2.5:7b"
```

---

## 支持的 AI 后端

`api_base` 是唯一的后端切换机制。修改该字段即可在不同后端之间切换，无需改动其他配置。

| 后端 | api_base | 适用场景 |
|------|----------|----------|
| **本地 vLLM** | `http://<局域网IP>:8000/v1` | 内网部署，数据不出本地，延迟低 |
| **DeepSeek API** | `https://api.deepseek.com/v1` | 云端 API，模型能力最强 |
| **Ollama** | `http://localhost:11434/v1` | 本地笔记本，轻量部署 |

管线通过 OpenAI SDK 统一调用，只要后端提供 OpenAI 兼容的 API 接口即可接入。

---

## AI 注释双模式

Step 05（`05_annotate_major.py`）采用双模式注释策略：

### 主模式：LLM 注释

1. 对每个聚类运行 Wilcoxon rank-sum 检验，提取 top N marker 基因
2. 使用 `build_annotation_prompt()` 组装提示词，包含组织、物种、每个聚类的 marker 基因列表
3. 调用 `ai_query()` 向 LLM 请求 JSON 格式的注释结果
4. LLM 返回结构：`{"cluster_id": {"cell_type": "...", "state": "...", "subtype": "...", "confidence": "...", "reasoning": "..."}}`

对于聚类数量较多的情况，注释会分块提交，以防止单次 prompt 超出模型的长度限制。

### 备用模式：score_genes

当 LLM 注释失败时，自动回退到 `score_genes` 方法：对每个预定义的细胞类型 marker 集计算平均表达得分，将聚类分配给得分最高的细胞类型。这种方式不需要外部 API 调用，完全离线运行。

---

## 回退链（Fallback Chain）

这是管线鲁棒性的关键设计。Step 05 的回退流程如下：

```
尝试 LLM 注释
    │
    ├─ ai_query() 抛出异常（网络错误/超时/API 拒绝）
    │   └─→ 记录警告，回退到 score_genes
    │
    ├─ LLM 返回 None（某些 vLLM 部署偶发空响应）
    │   └─→ 自动重试最多 3 次，仍为 None 则回退
    │
    ├─ 响应不是合法 JSON
    │   └─→ 记录警告 + 原始响应前 500 字符，回退
    │
    └─ 响应缺少必要字段（cell_type/state/subtype/confidence/reasoning）
        └─→ 记录哪个聚类缺了哪些字段，回退
```

所有分支均返回 `None`，调用方检测到 `None` 后调用 `score_genes_mode()` 完成注释。

**关键原则**：AI 注释从不阻塞管线。任何层面的失败都会触发优雅回退。

---

## AI 调用器（ai_caller.py）

`scripts/ai_caller.py` 是统一的 LLM 调用模块。

### ai_query() 接口

```python
def ai_query(system_prompt: str, user_prompt: str, cfg) -> str:
```

- 使用 OpenAI SDK 作为统一客户端
- `cfg` 为鸭子类型，只需提供 `api_base`、`model`、`api_key`、`max_tokens`、`temperature` 等属性
- 不导入 scanpy 或任何生物学分析库，保持调用层纯净

### 关键特性

- **思考模式**：通过 `thinking_enabled` / `reasoning_effort` 控制 DeepSeek 思考模式。开启时不传 `temperature`（API 忽略但语义更清晰）
- **自动重试**：某些 vLLM 部署会偶发返回 `content=None`，最多重试 3 次（指数退避：1s, 2s, 4s）
- **密钥回退**：`api_key` 为空时自动使用 `LLM_API_KEY` 环境变量。详见 [API Key 安全存储](#api-key-安全存储)

---

## Prompt 模板（ai_prompts.py）

`scripts/ai_prompts.py` 存放所有 LLM 提示词模板。

### 核心组件

| 组件 | 说明 |
|------|------|
| `ANNOTATION_SYSTEM_PROMPT` | 系统级角色提示词，定义生物学专家身份和 JSON 输出格式要求 |
| `ANNOTATION_USER_PROMPT_TEMPLATE` | 用户提示词模板，包含组织、物种、marker 基因占位符 |
| `build_annotation_prompt()` | 自动运行 marker 基因检测并组装完整提示词，返回 `(sys_prompt, user_prompt)` 二元组 |

### 提示词设计原则

- 提示词与调用逻辑分离（提示词在 `ai_prompts.py`，调用在 `ai_caller.py`）
- 要求 LLM 只返回 JSON，不包含 markdown 格式或额外解释
- 每条注释包含 `cell_type`、`state`、`subtype`、`confidence`、`reasoning` 五个字段
- `build_annotation_prompt()` 可自动计算或复用已有的 `rank_genes_groups` 结果

### 预留提示词

文件中预留了以下提示词桩（尚未实现）：

- `PARAM_SUGGEST_PROMPT`：根据数据特征建议 QC 参数
- `QC_REVIEW_PROMPT`：审查 QC 结果
- `DEG_DESIGN_PROMPT`：建议差异表达对比设计
- `INTERPRETATION_PROMPT`：解读差异表达或富集结果

---

## API Key 安全存储

API 密钥应避免硬编码在配置文件中，防止意外提交到 git。

### 推荐方式：`.env` 文件（本项目支持）

1. 复制模板文件：
   ```bash
   cp .env.example .env
   ```
2. 编辑 `.env`，填入从 [DeepSeek Platform](https://platform.deepseek.com/api_keys) 获取的 key：
   ```bash
   LLM_API_KEY=sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   ```
3. 确保 `.env` 已被 `.gitignore` 排除（默认已添加）。
4. `config.py` 中保持 `api_key = ""`（默认值），`ai_caller` 会自动读取 `LLM_API_KEY` 环境变量。

### 备选方式：直接设置环境变量

```bash
# 临时（当前终端）
export LLM_API_KEY=sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# 持久化（写入 shell 配置）
echo 'export LLM_API_KEY=sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' >> ~/.bashrc
source ~/.bashrc
```

### 优先级（高 → 低）

1. `config.py` 中 `cfg.ai.api_key`（显式指定）
2. 环境变量 `LLM_API_KEY`
3. 兜底值 `"not-needed"`（仅用于本地测试 vLLM 等无需认证的后端）

---

## 涉及的文件

| 文件 | 路径 | 作用 |
|------|------|------|
| AI 配置 | `config.py` → `AIConfig` | 控制 AI 功能开关和后端选择 |
| AI 调用器 | `scripts/ai_caller.py` | 封装 OpenAI SDK，提供统一 `ai_query()` 接口 |
| Prompt 模板 | `scripts/ai_prompts.py` | 存放所有提示词模板和构建函数 |
| 注释步骤 | `scripts/05_annotate_major.py` | 双模式注释的执行入口，包含完整回退链 |
| 亚型步骤 | `scripts/06_subcluster.py` | 调用 `ai_query()` 对亚聚类进行再注释 |
