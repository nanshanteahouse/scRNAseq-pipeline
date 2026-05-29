#!/usr/bin/env python3
"""
utils.py — 管线通用工具函数
==============================

集中存放各步骤共享的工具函数:
  - safe_write: WSL 兼容的 h5ad 安全保存
  - safe_plot: 容错的 matplotlib 绘图包装
  - setup_logger: 统一日志配置
"""

import os
import sys
import shutil
import logging
from pathlib import Path
from typing import Optional

import scanpy as sc


def safe_write(adata: sc.AnnData, target: str,
               tmpdir: str = "/tmp/scRNAseq_pipeline",
               compression: str = "gzip") -> None:
    """
    安全写入 h5ad 文件，避免 WSL /mnt 挂载的文件锁定问题。

    策略: 先写入 /tmp，再 mv 到目标路径。
    mv 是原子操作（在同一文件系统内），确保不会留下损坏的中间文件。

    GSE169109 项目的 safe_write 模式在此复用。

    参数:
        adata: AnnData 对象
        target: 目标 .h5ad 路径
        tmpdir: 临时目录
        compression: h5py 压缩方式 ('gzip' | 'lzf' | None)
    """
    os.makedirs(tmpdir, exist_ok=True)
    tmp_path = os.path.join(tmpdir, os.path.basename(target))
    adata.write(tmp_path, compression=compression)
    shutil.move(tmp_path, target)

    size_mb = os.path.getsize(target) / 1e6
    logger = logging.getLogger(__name__)
    logger.info("Saved %s (%.1f MB)", os.path.basename(target), size_mb)


def safe_plot(func, *args, **kwargs):
    """
    容错的 scanpy 绘图包装。

    某些 scanpy 绘图函数在某些版本组合下可能因 matplotlib 兼容性崩溃。
    本函数捕获异常并记录警告，避免整个步骤因此中断。

    用法:
        safe_plot(sc.pl.umap, adata, color='stage', show=False, save='_stage.pdf')
    """
    logger = logging.getLogger(__name__)
    try:
        return func(*args, **kwargs)
    except Exception as e:
        logger.warning("Plot failed (skipped): %s", e)
        return None


def setup_logger(name: str, log_file: str,
                 level: int = logging.INFO) -> logging.Logger:
    """
    统一配置日志: 同时输出到 stdout 和文件。

    格式模仿 GSE169109 项目:
        14:30:00 | INFO    | 消息内容

    参数:
        name: logger 名称
        log_file: 日志文件路径
        level: 日志级别

    返回:
        配置好的 logger 实例
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # 防止重复添加 handler
    if logger.handlers:
        return logger

    formatter = logging.Formatter(
        '%(asctime)s | %(levelname)-7s | %(message)s',
        datefmt='%H:%M:%S',
    )

    # stdout handler
    stdout_h = logging.StreamHandler(sys.stdout)
    stdout_h.setFormatter(formatter)
    logger.addHandler(stdout_h)

    # file handler
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    file_h = logging.FileHandler(log_file, mode='w')
    file_h.setFormatter(formatter)
    logger.addHandler(file_h)

    return logger


def resolve_config(config_path: Optional[str] = None):
    """
    解析 --config 参数，返回配置模块的 CFG 对象。

    所有步骤脚本统一使用本函数加载配置。
    """
    if config_path is None:
        # 默认寻找父目录的 config.py
        config_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "config.py",
        )

    config_path = os.path.abspath(config_path)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"配置文件不存在: {config_path}")

    import importlib.util
    spec = importlib.util.spec_from_file_location("pipeline_config", config_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pipeline_config"] = mod
    spec.loader.exec_module(mod)
    return mod.CFG
