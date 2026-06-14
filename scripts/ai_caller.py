#!/usr/bin/env python3
"""
ai_caller.py — 统一 LLM 调用模块
==================================

集中管理所有 AI API 调用，提供统一的接口模式。
所有步骤脚本通过本模块调用 LLM，避免重复的客户端创建和错误处理。

设计原则:
  - 使用 OpenAI SDK 作为统一客户端（兼容 OpenAI / Azure / 兼容 API）
  - cfg 为鸭子类型，只需提供 api_base, model, api_key, max_tokens, temperature 属性
  - 不在此模块中导入 scanpy 或任何生物学分析库

用法:
    from ai_caller import ai_query
    from ai_prompts import ANNOTATION_SYSTEM_PROMPT

    result = ai_query(ANNOTATION_SYSTEM_PROMPT, "User query...", cfg)
"""

import os
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def ai_query(system_prompt: str, user_prompt: str, cfg) -> str:
    """
    统一的 LLM 查询接口。

    使用 OpenAI SDK 构造聊天补全请求。cfg 是任意提供以下属性的对象:
      - api_base (str):      API 端点 URL
      - model (str):         模型名称 (如 'gpt-4o', 'claude-3-opus-20240229')
      - api_key (str):       API 密钥（可空，由 SDK 自动回退环境变量）
      - max_tokens (int):    最大输出 token 数
      - temperature (float): 采样温度 (0.0 ~ 2.0)

    参数:
        system_prompt: 系统角色提示词
        user_prompt:   用户输入/任务提示词
        cfg:           配置对象（鸭子类型，只需包含上述属性）

    返回:
        模型生成的文本内容

    示例:
        >>> resp = ai_query("你是一个生物学专家", "注释以下聚类...", llm_cfg)
        >>> print(resp)
        '{"0": {"cell_type": "T cell", ...}}'
    """
    from openai import OpenAI

    client = OpenAI(
        api_key=cfg.api_key or os.getenv("LLM_API_KEY", "not-needed"),
        base_url=cfg.api_base,
    )

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user",   "content": user_prompt},
    ]

    # Compute effective max_tokens: if reasoning_budget is set, add it to max_tokens
    reasoning_budget = getattr(cfg, "reasoning_budget", 0) or 0
    effective_max_tokens = cfg.max_tokens + reasoning_budget

    # Retry loop: some vLLM deployments return content=None transiently
    max_retries = 3
    for attempt in range(max_retries):
        resp = client.chat.completions.create(
            model=cfg.model,
            messages=messages,
            max_tokens=effective_max_tokens,
            temperature=cfg.temperature,
            timeout=getattr(cfg, "timeout", None),
        )
        content = resp.choices[0].message.content
        if content is not None:
            return content
        if attempt < max_retries - 1:
            wait = 2 ** attempt
            print(f"[ai_caller] Empty response, retrying in {wait}s...", file=sys.stderr)
            time.sleep(wait)

    return None
