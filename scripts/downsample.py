#!/usr/bin/env python3
"""
Downsampling — 大型 scRNA-seq 数据集的细胞降采样
===================================================
当数据集过大导致下游步骤 OOM 时，在任意 h5ad checkpoint 之间插入此脚本，
减少细胞数以控制内存使用。

支持三种降采样策略:
  1. stratified  (默认): 按样本分层采样，保持各样本比例 → 适合有 sample 列的数据
  2. random:             完全随机采样 → 适合无分组的简单降采样
  3. max_per_sample:     每个样本最多保留 N 个细胞 → 适合样本大小极不均衡的数据

使用方法:
  # 插入在 pipeline 步骤之间:
  python run_pipeline.py --step 0 --config config_large.py
  ./venv/bin/python scripts/downsample.py --config config_large.py --target-total 50000
  python run_pipeline.py --step 1 --config config_large.py

  # 指定输入输出 checkpoint:
  ./venv/bin/python scripts/downsample.py \\
      --input results/h5ad/00_raw.h5ad \\
      --output results/h5ad/00_raw_downsampled.h5ad \\
      --strategy stratified --target-total 30000

  # 超大文件使用 backed 模式（低内存读取）:
  ./venv/bin/python scripts/downsample.py \\
      --config config_large.py --backed --target-total 50000

  # 每个样本最多 3000 细胞:
  ./venv/bin/python scripts/downsample.py \\
      --config config_large.py --strategy max_per_sample --max-per-sample 3000
"""
import sys, os, time, argparse
import numpy as np
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import setup_logger, resolve_config, safe_write

# 与 run_pipeline.py CHECKPOINT_FILES 同步，副本避免跨目录导入
CHECKPOINT_FILES = [
    "00_raw.h5ad", "01_doublet.h5ad", "02_qc.h5ad",
    "03_integrated.h5ad", "04_clustered.h5ad", "05_annotated.h5ad",
    "05_annotated.h5ad", "05_annotated.h5ad", "04_clustered.h5ad",
    "marker_genes_per_group.csv", "05_annotated.h5ad",
]


def _check_sample_col(adata: sc.AnnData, sample_key: str, log) -> str:
    """查找可用的样本分组列。返回实际使用的列名或 None。"""
    if sample_key and sample_key in adata.obs:
        return sample_key
    # 尝试常见列名
    for candidate in ['sample', 'Sample', 'samples', 'batch', 'Batch', 'stage', 'Stage']:
        if candidate in adata.obs:
            log.info("使用 '%s' 作为分组列（未找到 '%s'）", candidate, sample_key)
            return candidate
    return None


def downsample_random(adata: sc.AnnData, target: int, rng: np.random.RandomState,
                      log) -> sc.AnnData:
    """完全随机采样 target 个细胞。"""
    n_cells = adata.n_obs
    if target >= n_cells:
        log.info("target_total (%d) >= 当前细胞数 (%d)，无需降采样", target, n_cells)
        return adata
    idx = rng.choice(n_cells, size=target, replace=False)
    idx.sort()
    log.info("随机采样: %d → %d 细胞 (%.1f%%)", n_cells, target, 100 * target / n_cells)
    return adata[idx].copy()


def downsample_stratified(adata: sc.AnnData, target: int, sample_key: str,
                          rng: np.random.RandomState, log) -> sc.AnnData:
    """按样本分层采样，保持各样本比例。"""
    n_cells = adata.n_obs
    if target >= n_cells:
        log.info("target_total (%d) >= 当前细胞数 (%d)，无需降采样", target, n_cells)
        return adata

    counts = adata.obs[sample_key].value_counts()
    log.info("分层采样, 分组=%s, 目标总细胞=%d", sample_key, target)
    for s, c in counts.items():
        log.info("  样本 %s: %d 细胞 (%.1f%%)", s, c, 100 * c / n_cells)

    # 按比例分配 target
    fractions = counts / n_cells
    per_sample_targets = (fractions * target).astype(int)
    # 处理余数 — 从余数最大的样本补 1
    remainder = target - per_sample_targets.sum()
    if remainder > 0:
        sorted_idx = np.argsort((fractions * target) - per_sample_targets)[::-1]
        for i in range(remainder):
            per_sample_targets.iloc[int(sorted_idx[i])] += 1

    # 每个样本分别采样
    indices = []
    for sample_name in counts.index:
        mask = adata.obs[sample_key] == sample_name
        sample_idx = np.where(mask)[0]
        n_sample = len(sample_idx)
        t = min(per_sample_targets[sample_name], n_sample)
        if t < n_sample:
            chosen = rng.choice(sample_idx, size=t, replace=False)
        else:
            chosen = sample_idx
        indices.append(chosen)

    idx = np.concatenate(indices)
    idx.sort()
    log.info("分层采样: %d → %d 细胞 (%.1f%%)", n_cells, len(idx), 100 * len(idx) / n_cells)
    return adata[idx].copy()


def downsample_max_per_sample(adata: sc.AnnData, max_per: int, sample_key: str,
                              rng: np.random.RandomState, log) -> sc.AnnData:
    """每个样本最多保留 max_per 个细胞。"""
    counts = adata.obs[sample_key].value_counts()
    log.info("按样本封顶, 每样本最多 %d 细胞", max_per)

    indices = []
    for sample_name in counts.index:
        mask = adata.obs[sample_key] == sample_name
        sample_idx = np.where(mask)[0]
        n_sample = len(sample_idx)
        if n_sample > max_per:
            chosen = rng.choice(sample_idx, size=max_per, replace=False)
            log.info("  样本 %s: %d → %d (截断 %d)", sample_name, n_sample, max_per, n_sample - max_per)
        else:
            chosen = sample_idx
            log.info("  样本 %s: %d (不变)", sample_name, n_sample)
        indices.append(chosen)

    idx = np.concatenate(indices)
    idx.sort()
    log.info("封顶采样: %d → %d 细胞 (%.1f%%)", adata.n_obs, len(idx), 100 * len(idx) / adata.n_obs)
    return adata[idx].copy()


def estimate_memory_gb(adata: sc.AnnData) -> float:
    """粗略估计 AnnData 在内存中的大小 (GB)。"""
    total = 0.0
    # X matrix
    if hasattr(adata, 'X') and adata.X is not None:
        if sp.issparse(adata.X):
            # CSR: data + indices + indptr
            total += adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes
        else:
            total += adata.X.nbytes
    # obs
    for col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        if dtype == object:
            continue  # string columns are harder to estimate
        total += adata.obs[col].values.nbytes if hasattr(adata.obs[col].values, 'nbytes') else 0
    # var
    for col in adata.var.columns:
        dtype = adata.var[col].dtype
        if dtype == object:
            continue
        total += adata.var[col].values.nbytes if hasattr(adata.var[col].values, 'nbytes') else 0
    # layers
    if hasattr(adata, 'layers'):
        for layer_name in adata.layers.keys():
            layer = adata.layers[layer_name]
            if sp.issparse(layer):
                total += layer.data.nbytes + layer.indices.nbytes + layer.indptr.nbytes
            elif layer is not None:
                total += layer.nbytes
    # obsm
    if hasattr(adata, 'obsm'):
        for key in adata.obsm.keys():
            arr = adata.obsm[key]
            if hasattr(arr, 'nbytes'):
                total += arr.nbytes
    # varm
    if hasattr(adata, 'varm'):
        for key in adata.varm.keys():
            arr = adata.varm[key]
            if hasattr(arr, 'nbytes'):
                total += arr.nbytes
    # uns (skip, too heterogeneous)
    return total / (1024 ** 3)


def main():
    t0 = time.time()

    parser = argparse.ArgumentParser(
        description="scRNA-seq 细胞降采样 — 配合 pipeline 使用",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--config", default=None,
                        help="配置文件路径（提供后自动读取 h5ad_dir 等路径）")
    parser.add_argument("--input", "--input-h5ad", dest="input_h5ad", default=None,
                        help="输入 h5ad 路径（默认: config 中指定步骤的 checkpoint）")
    parser.add_argument("--output", "--output-h5ad", dest="output_h5ad", default=None,
                        help="输出 h5ad 路径（默认: 覆盖输入文件）")
    parser.add_argument("--step", type=int, default=0,
                        help="从哪一步的 checkpoint 读取（默认 0 = 00_raw.h5ad）")

    strategy_group = parser.add_mutually_exclusive_group()
    strategy_group.add_argument("--strategy", default="stratified",
                                choices=["random", "stratified", "max_per_sample"],
                                help="降采样策略（默认: stratified）")
    strategy_group.add_argument("--random", action="store_true",
                                help="等价于 --strategy random")

    target_group = parser.add_mutually_exclusive_group()
    target_group.add_argument("--target-total", type=int, default=None,
                              help="目标总细胞数（管道模式从 config.downsample_target 读取）")
    target_group.add_argument("--target-fraction", type=float, default=None,
                              help="保留比例 (0.0-1.0)")
    target_group.add_argument("--max-per-sample", type=int, default=None,
                              help="每个样本最多保留 N 细胞（仅 max_per_sample 策略）")

    parser.add_argument("--sample-key", type=str, default="sample",
                        help="obs 中样本名列（默认: sample）")
    parser.add_argument("--random-seed", type=int, default=42,
                        help="随机种子（默认: 42）")
    parser.add_argument("--backed", action="store_true",
                        help="使用 backed='r' 模式读取 h5ad（低内存，适用于超大数据集）")
    parser.add_argument("--dry-run", action="store_true",
                        help="仅预览统计信息，不执行降采样")
    parser.add_argument("--overwrite", action="store_true",
                        help="直接覆盖输入文件（而非写入新文件）")
    parser.add_argument("--to-config", action="store_true",
                        help="将降采样结果作为后续步骤的输入（更新 config checkpoint 引用）")

    args = parser.parse_args()

    # ── 解析策略别名 ──
    strategy = args.strategy
    if args.random:
        strategy = "random"
    if strategy == "max_per_sample":
        if args.max_per_sample is None:
            parser.error("--strategy max_per_sample 需要 --max-per-sample 参数")
        target_per_sample = args.max_per_sample
    elif args.max_per_sample is not None:
        parser.error("--max-per-sample 仅用于 --strategy max_per_sample")

    # ── 加载配置 ──
    CFG = None
    if args.config:
        CFG = resolve_config(args.config)
        CFG.resolve_paths()

    # ── 确定输入/输出路径 ──
    if args.input_h5ad:
        input_path = args.input_h5ad
    elif CFG is not None:
        step_idx = min(args.step, len(CHECKPOINT_FILES) - 1)
        input_path = os.path.join(CFG.h5ad_dir, CHECKPOINT_FILES[step_idx])
    else:
        parser.error("请提供 --config 或 --input")

    if args.output_h5ad:
        output_path = args.output_h5ad
    elif args.overwrite or args.config is None:
        output_path = input_path
    elif CFG is not None:
        # 默认生成 *_downsampled.h5ad 同级文件
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}_downsampled{ext}"
    else:
        output_path = input_path

    # ── 从 config 读取降采样参数（仅当未通过 CLI 指定时） ──
    if CFG is not None:
        if args.target_total is None and args.target_fraction is None and args.max_per_sample is None:
            if CFG.downsample_target is not None:
                args.target_total = CFG.downsample_target
                if CFG.downsample_strategy:
                    strategy = CFG.downsample_strategy
                if strategy == "max_per_sample" and CFG.downsample_max_per_sample is not None:
                    args.max_per_sample = CFG.downsample_max_per_sample
                args.random_seed = CFG.downsample_random_seed
            else:
                # 配置未启用降采样 — 跳过（管道模式无操作）
                print("[downsample] 跳过: downsample_target 未配置")
                return

    # ── 设置日志 ──
    if CFG is not None:
        log_dir = CFG.log_dir
    else:
        # 无 config 时，在输出文件同目录下写日志
        log_dir = os.path.join(os.path.dirname(os.path.abspath(output_path)), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    log = setup_logger("downsample", os.path.join(log_dir, "downsample.log"))
    log.info("=" * 60)
    log.info("降采样 — Downsampling")
    log.info("=" * 60)
    log.info("输入: %s", input_path)
    log.info("输出: %s", output_path)
    log.info("策略: %s", strategy)
    if args.target_total:
        log.info("目标总细胞: %d", args.target_total)
    elif args.target_fraction:
        log.info("保留比例: %.2f", args.target_fraction)
    elif strategy == "max_per_sample":
        log.info("每样本上限: %d", args.max_per_sample)

    if not os.path.exists(input_path):
        log.error("输入文件不存在: %s", input_path)
        sys.exit(1)
    if args.dry_run:
        log.info("[DRY RUN] 仅预览，不执行降采样")

    # ── 读取数据 ──
    log.info("读取数据...")
    read_t = time.time()
    try:
        adata = sc.read(input_path, backed='r' if args.backed else None)
    except Exception as e:
        log.error("读取失败: %s", e)
        sys.exit(1)
    log.info("  细胞: %d × 基因: %d", adata.n_obs, adata.n_vars)
    log.info("  读取耗时: %.1fs", time.time() - read_t)

    # 如果 backed，加载到内存
    if args.backed:
        log.info("  backed 模式 — 转换为内存中 AnnData...")
        adata = adata.to_memory()

    # 预估内存
    est_gb = estimate_memory_gb(adata)
    log.info("  预估内存: %.2f GB (仅此对象)", est_gb)

    # ── 确定采样目标 ──
    n_cells = adata.n_obs
    if args.target_fraction is not None:
        target_total = int(n_cells * args.target_fraction)
        log.info("保留比例 %.2f → 目标 %d 细胞", args.target_fraction, target_total)
    elif args.target_total is not None:
        target_total = args.target_total
    else:
        target_total = n_cells  # max_per_sample handles itself

    # ── 预览统计 ──
    if strategy != "max_per_sample":
        log.info("当前: %d 细胞 → 目标: %d 细胞 (减少 %.1f%%)",
                  n_cells, target_total,
                  (1 - target_total / n_cells) * 100 if target_total < n_cells else 0)
    else:
        log.info("当前: %d 细胞 → 每样本上限: %d", n_cells, args.max_per_sample)

    sample_col = _check_sample_col(adata, args.sample_key, log)
    if sample_col:
        log.info("各分组细胞数:")
        counts = adata.obs[sample_col].value_counts()
        for s, c in counts.items():
            log.info("  %-30s %6d cells", s, c)

    if args.dry_run:
        log.info("[DRY RUN] 预览完成。移除 --dry-run 来执行降采样。")
        log.info("提示: 如果预估内存不足，建议 target_total ≤ %d",
                 int(n_cells * 0.3))
        return

    # ── 执行降采样 ──
    rng = np.random.RandomState(args.random_seed)
    log.info("执行降采样 (策略=%s)...", strategy)

    if strategy == "random":
        adata = downsample_random(adata, target_total, rng, log)
    elif strategy == "stratified":
        if sample_col is None:
            log.warning("未找到分组列，回退到随机采样")
            adata = downsample_random(adata, target_total, rng, log)
        else:
            adata = downsample_stratified(adata, target_total, sample_col, rng, log)
    elif strategy == "max_per_sample":
        if sample_col is None:
            log.warning("未找到分组列，回退到随机采样")
            adata = downsample_random(adata, target_total, rng, log)
        else:
            adata = downsample_max_per_sample(adata, args.max_per_sample, sample_col, rng, log)

    # 可选 float32 节省内存
    if CFG and getattr(CFG, 'use_float32', False):
        if sp.issparse(adata.X):
            adata.X = adata.X.astype('float32', copy=False)
        else:
            adata.X = adata.X.astype('float32')
        log.info("X 精度已转换为 float32")

    # ── 写入 ──
    log.info("保存到: %s", output_path)
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    if output_path == input_path and CFG:
        # 覆盖原始 checkpoint — 使用 safe_write
        safe_write(adata, output_path)
    else:
        # 新文件 — safe_write 是安全的
        safe_write(adata, output_path)

    # ── 摘要 ──
    elapsed = time.time() - t0
    new_est = estimate_memory_gb(adata)
    log.info("=" * 60)
    log.info("降采样完成!")
    log.info("  %d → %d 细胞 (%d 基因)", n_cells, adata.n_obs, adata.n_vars)
    log.info("  预估内存: %.2f GB → %.2f GB", est_gb, new_est)
    log.info("  耗时: %.1fs (%.1fmin)", elapsed, elapsed / 60)
    log.info("  输出: %s", output_path)
    log.info("=" * 60)

    if CFG and args.to_config:
        log.info("提示: 将降采样后的文件重命名为原始 checkpoint 名可无缝衔接后续步骤")
        log.info("  mv %s %s", output_path, input_path)

    print(f"\n✅ 降采样完成: {n_cells} → {adata.n_obs} 细胞")
    print(f"   输出: {output_path}")
    print(f"   耗时: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
