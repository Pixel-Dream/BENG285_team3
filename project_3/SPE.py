#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_sigprofiler.py

用 Python API 调用 SigProfilerExtractor 完整跑批：
  • Poisson 重抽样  
  • 多次 NMF  
  • 聚类稳定性分析  
  • 自动选 K  
依赖：pip install sigprofilerextractor
"""

import os
import sys
from SigProfilerExtractor import sigpro as sp

def main():
    # --------------- 配置区 ---------------
    # 输入矩阵文件（matrix 格式）
    input_matrix = "SBS/TCGA.SBS96.all"
    # 输出目录
    output_dir   = "SPx_output"
    # signature 数目搜索范围
    min_k = 1
    max_k = 25
    # NMF 复刻次数
    nmf_replicates = 100
    # 并行 CPU 核数
    cpu = 8
    # 通道类型（对应文件里的类别数，比如 SBS96、SBS384…）
    context = "SBS96"
    # ---------------------------------------

    # 先检查输入文件
    if not os.path.exists(input_matrix):
        print(f"❌ 找不到输入矩阵：{input_matrix}", file=sys.stderr)
        sys.exit(1)

    print("▶ 开始运行 SigProfilerExtractor …")
    sp.sigProfilerExtractor(
        input_type="matrix",            # 输入类型：matrix / vcf / maf
        output=output_dir,              # 输出目录
        input_data=input_matrix,        # 矩阵文件路径
        reference_genome=context,       # 通道上下文
        minimum_signatures=min_k,       # K 最小值
        maximum_signatures=max_k,       # K 最大值
        nmf_replicates=nmf_replicates,  # 每个 K 的 NMF 重复次数
        cpu=cpu,                        # 并行线程数
        # 以下为可选高阶参数，按需取消注释：
        # minimum_nmf_iterations=1000,
        # maximum_nmf_iterations=200000,
        # nmf_test_convergence=1000,
        # nmf_tolerance=1e-08,
    )
    print("✅ SigProfilerExtractor 运行完毕。结果保存在：", output_dir)

if __name__ == "__main__":
    main()
