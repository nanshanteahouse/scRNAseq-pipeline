#!/usr/bin/env python3
"""
run_pipeline.py — 通用 scRNA-seq 管线主控
==========================================

功能:
  - 按顺序执行全部 10 步分析
  - 指定单步或步骤范围运行
  - 从断点恢复（检查 h5ad checkpoint 是否存在）
  - 列出所有可用步骤

用法:
    python run_pipeline.py                      # 全部顺序执行
    python run_pipeline.py --step 3             # 只跑第 3 步
    python run_pipeline.py --steps 3-6          # 跑 3~6 步
    python run_pipeline.py --steps 1,3,5        # 跑第 1, 3, 5 步
    python run_pipeline.py --resume             # 从第一个未完成的步骤继续
    python run_pipeline.py --list               # 列出所有步骤
    python run_pipeline.py --config my_config.py # 使用自定义配置

checkpoint 依赖链:
    00_raw.h5ad  ← 00_load.py
    01_qc.h5ad   ← 01_qc.py
    02_normalized.h5ad ← 02_normalize_hvg.py
    03_harmony.h5ad    ← 03_pca_harmony.py
    04_clustered.h5ad  ← 04_cluster_umap.py (步骤 05-09 读取此 checkpoint)
    05_final.h5ad      ← 08_trajectory.py
"""

import sys
import os
import subprocess
import argparse


# ── 步骤注册表 ──────────────────────────────────────────────────────────
# 每步: (序号, 脚本名, 描述)
STEPS = [
    ("00", "00_load.py",          "加载原始数据 → 00_raw.h5ad"),
    ("01", "01_qc.py",            "QC 指标 + Scrublet + 过滤 → 01_qc.h5ad"),
    ("02", "02_normalize_hvg.py", "归一化 + HVG 选择 → 02_normalized.h5ad"),
    ("03", "03_pca_harmony.py",   "PCA + Harmony 批次校正 → 03_harmony.h5ad"),
    ("04", "04_cluster_umap.py",  "邻居图 + UMAP + 多分辨率 Leiden → 04_clustered.h5ad"),
    ("05", "05_annotate.py",      "细胞类型自动注释 + 子聚类"),
    ("06", "06_exploratory.py",   "细胞组成 + QC + 标记基因探索"),
    ("07", "07_markers_de.py",    "标记基因 + 组间差异表达"),
    ("08", "08_trajectory.py",    "PAGA + DPT + 分支分析 → 05_final.h5ad"),
    ("09", "09_enrichment.py",    "GO/KEGG 通路富集分析"),
]

# 每步对应的 checkpoint 文件名
# 步骤 05-08 读取 04_clustered.h5ad，不写新的 checkpoint
# 步骤 09 读取 marker_genes_per_group.csv (Step 07 产出), 不写 h5ad
CHECKPOINT_FILES = [
    "00_raw.h5ad",
    "01_qc.h5ad",
    "02_normalized.h5ad",
    "03_harmony.h5ad",
    "04_clustered.h5ad",
    "04_clustered.h5ad",   # step 05 reads
    "04_clustered.h5ad",   # step 06 reads
    "04_clustered.h5ad",   # step 07 reads
    "04_clustered.h5ad",   # step 08 reads
    "04_clustered.h5ad",   # step 09 reads (no h5ad output)
]

# 哪些步骤会输出新的 checkpoint（用于 --resume 判断）
STEPS_WRITE_CHECKPOINT = {0, 1, 2, 3, 4, 8}


def find_first_incomplete(h5ad_dir: str) -> int:
    """
    扫描 checkpoint 目录，找到第一个未完成的步骤。

    返回:
        第一个未完成步骤的索引，如果全部完成则返回 len(STEPS)
    """
    for i in range(len(STEPS)):
        if i not in STEPS_WRITE_CHECKPOINT:
            continue  # 跳过不写 checkpoint 的步骤
        ckpt = os.path.join(h5ad_dir, CHECKPOINT_FILES[i])
        if not os.path.exists(ckpt):
            return i
    return len(STEPS)


def parse_step_range(spec: str) -> list:
    """
    解析步骤范围字符串。

    支持格式:
      - "3-6"  → [3, 4, 5, 6]
      - "1,3,5" → [1, 3, 5]
      - "3"    → [3]
    """
    if "-" in spec:
        a, b = map(int, spec.split("-"))
        return list(range(a, b + 1))
    else:
        return [int(s) for s in spec.split(",")]


def main():
    parser = argparse.ArgumentParser(
        description="通用 scRNA-seq 分析管线主控",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--steps", type=str,
                       help="步骤范围 (如 3-6) 或列表 (如 1,3,5)")
    group.add_argument("--step", type=int,
                       help="运行单个步骤（1-based）")
    group.add_argument("--resume", action="store_true",
                       help="从第一个未完成的 checkpoint 继续")
    parser.add_argument("--list", action="store_true",
                        help="列出所有步骤")
    parser.add_argument("--config", type=str, default="config.py",
                        help="配置文件路径（默认: config.py）")
    args = parser.parse_args()

    # ── --list 模式 ────────────────────────────────────────────────
    if args.list:
        print("通用 scRNA-seq 分析管线 — 步骤列表")
        print("=" * 60)
        for num, script, desc in STEPS:
            ckpt = CHECKPOINT_FILES[STEPS.index((num, script, desc))]
            print(f"  [{num}] {desc}")
            print(f"        脚本: {script}  |  checkpoint: {ckpt}")
        print(f"\n用法: python {os.path.basename(__file__)} --step N")
        return

    # ── 加载配置 ───────────────────────────────────────────────────
    # 动态导入用户指定的配置文件
    config_path = os.path.abspath(args.config)
    if not os.path.exists(config_path):
        print(f"[run] 错误: 配置文件不存在: {config_path}")
        sys.exit(1)

    import importlib.util
    spec = importlib.util.spec_from_file_location("pipeline_config", config_path)
    cfg_module = importlib.util.module_from_spec(spec)
    sys.modules["pipeline_config"] = cfg_module
    spec.loader.exec_module(cfg_module)
    CFG = cfg_module.CFG
    CFG.resolve_paths()

    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    python_exe = sys.executable

    # ── 解析步骤范围 ──────────────────────────────────────────────
    if args.resume:
        start = find_first_incomplete(CFG.h5ad_dir)
        if start >= len(STEPS):
            print("[run] ✅ 所有步骤已完成。")
            return
        step_indices = list(range(start, len(STEPS)))
        print(f"[run] 从步骤 [{STEPS[start][0]}] 恢复执行")
    elif args.steps:
        step_indices_raw = parse_step_range(args.steps)
        step_indices = [i - 1 for i in step_indices_raw]
        for i in step_indices:
            if i < 0 or i >= len(STEPS):
                print(f"[run] 错误: 无效步骤号 {i + 1}（有效范围: 1-{len(STEPS)}）")
                sys.exit(1)
    elif args.step:
        if args.step < 1 or args.step > len(STEPS):
            print(f"[run] 错误: 步骤号 {args.step} 超出范围 (1-{len(STEPS)})")
            sys.exit(1)
        step_indices = [args.step - 1]
    else:
        step_indices = list(range(len(STEPS)))

    # ── 执行步骤 ──────────────────────────────────────────────────
    for i in step_indices:
        num, script, desc = STEPS[i]
        script_path = os.path.join(scripts_dir, "scripts", script)

        if not os.path.exists(script_path):
            print(f"[run] 错误: 脚本不存在: {script_path}")
            sys.exit(1)

        print(f"\n{'=' * 60}")
        print(f"[run] 步骤 [{num}]: {desc}")
        print(f"{'=' * 60}")

        result = subprocess.run(
            [python_exe, script_path, f"--config={config_path}"],
            cwd=scripts_dir,
        )

        if result.returncode != 0:
            print(f"\n[run] ❌ 步骤 [{num}] 失败 (exit code={result.returncode})")
            print(f"[run] 修复问题后可通过以下命令继续:")
            print(f"      python {__file__} --resume --config {args.config}")
            sys.exit(1)

        print(f"[run] ✅ 步骤 [{num}] 完成。")

    print(f"\n{'=' * 60}")
    print(f"[run] 🎯 管线执行完毕。")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
