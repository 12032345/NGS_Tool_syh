import subprocess
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime
import argparse

def get_counts(bam_file, chrom, start, end, bin_size):
    """使用samtools计算指定区间内每个bin的符合条件的reads数"""
    bins = range(start, end, bin_size)
    counts = []
    
    # 检查samtools是否可用
    try:
        subprocess.run(['samtools', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("[ERROR] samtools not found in PATH")
        return list(bins), [0] * len(bins)
    
    for b_start in bins:
        b_end = b_start + bin_size
        bin_count = 0
        
        try:
            # 使用samtools view获取指定区域的reads
            cmd = [
                'samtools', 'view', '-c', 
                bam_file, 
                f'{chrom}:{b_start}-{b_end-1}'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            total_reads = int(result.stdout.strip())
            
            # 对于简化的计算，我们使用总reads数
            # 在实际应用中，更精确的过滤需要额外的samtools命令
            bin_count = total_reads
            
        except subprocess.CalledProcessError as e:
            # 如果区域不存在或samtools出错，计数为0
            bin_count = 0
        except ValueError:
            bin_count = 0
            
        counts.append(bin_count)
        
    return list(bins), counts

def plot_data(bins, counts, chrom, bin_size, title, filename, target_pos=None):
    """绘图并保存

    使用细竖线表示每个 bin（按光谱配色）。如果提供 `target_pos`，会在该位置画一条竖直虚线并标注原始坐标值。
    """
    import matplotlib.cm as cm
    from matplotlib.colors import Normalize

    plt.figure(figsize=(12, 5))
    # 使用 bin 中心作为每个竖线的位置，单位 Mb
    x_centers = (np.array(bins) + bin_size / 2.0) / 1e6

    counts_arr = np.array(counts)

    # 颜色映射：基于丰度（counts）映射颜色，低值偏蓝，高值偏红
    cmap = cm.get_cmap('RdYlBu_r')
    max_count = counts_arr.max() if counts_arr.size else 1
    norm = Normalize(vmin=0, vmax=max(max_count, 1))

    # 画细竖线（linewidth 很小，看起来像细柱），颜色由高度决定
    for x, h in zip(x_centers, counts_arr):
        if h > 0:
            c = cmap(norm(h))
            plt.vlines(x, 0, h, color=c, linewidth=0.9)
        else:
            # 对于 0 值使用浅灰色做占位
            plt.vlines(x, 0, 0, color=(0.9, 0.9, 0.9), linewidth=0.4)

    plt.title(title, fontsize=14)
    plt.xlabel(f"Chromosome {chrom} Position (Mb)", fontsize=12)
    plt.ylabel("Read Counts (Filtered)", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.3)

    # 如果给定目标位置，则绘制垂直虚线并标注原始坐标
    if target_pos is not None:
        x_target_mb = target_pos / 1e6
        # 竖线与标签半透明（alpha=0.6）
        plt.axvline(x=x_target_mb, color='red', linestyle='--', linewidth=1, alpha=0.6)
        # 在图顶端标注原始坐标值（不缩放到 Mb，显示整数坐标）
        ymax = counts_arr.max() if counts_arr.size else 1
        # 将文字放在竖线稍上方并倾斜90度以与竖线对齐
        plt.text(x_target_mb, ymax * 0.95, f"{int(target_pos)}", rotation=90,
                 va='top', ha='right', color=(1.0, 0.0, 0.0, 0.6), fontsize=10,
                 backgroundcolor=(1.0, 1.0, 1.0, 0.6))

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"✅ 成功生成图像: {filename}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='WGSmapping: WGS background and target enrichment plotting')
    # 适配bash脚本的参数调用方式
    parser.add_argument('--bam', required=True, help='Input BAM file path')
    parser.add_argument('--chromosome', required=True, help='Target chromosome (e.g., chr6)')
    parser.add_argument('--center', type=int, required=True, help='Center position (bp)')
    parser.add_argument('--step', type=int, default=100000, help='Step size for genome-wide analysis (bp)')
    parser.add_argument('--background', type=str, default='true', help='Perform background analysis (true/false)')
    parser.add_argument('--output', required=True, help='Output directory for plots')

    args = parser.parse_args()

    print("[INFO] === WORF-Seq 染色体比对分析 ===")
    print(f"[INFO] BAM文件: {args.bam}")
    print(f"[INFO] 染色体: {args.chromosome}")
    print(f"[INFO] 中心位置: {args.center}")
    print(f"[INFO] 步长: {args.step}")
    print(f"[INFO] 背景分析: {args.background}")

    # 1. 检查BAM文件
    bam_path = args.bam
    if not os.path.exists(bam_path):
        print(f"[ERROR] BAM文件不存在: {bam_path}")
        return
    
    if not os.path.exists(bam_path + ".bai"):
        print(f"[WARNING] 未找到索引文件 (.bai): {bam_path}.bai")

    # 检查BAM文件是否可读
    try:
        # 使用samtools view测试BAM文件
        result = subprocess.run(['samtools', 'view', '-c', bam_path], 
                              capture_output=True, text=True, check=True)
        total_reads = int(result.stdout.strip())
        print(f"[INFO] BAM文件包含 {total_reads:,} 条reads")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] 无法读取BAM文件: {e}")
        return
    except FileNotFoundError:
        print("[ERROR] samtools not found in PATH")
        return

    # 2. 输出目录设置
    out_dir = args.output
    os.makedirs(out_dir, exist_ok=True)
    
    # 获取文件基础名用于生成文件名
    bam_basename = os.path.splitext(os.path.basename(bam_path))[0]

    # 3. 设置目标染色体和位置
    target_chrom = args.chromosome
    target_pos = args.center
    wgs_bin = args.step
    
    # 转换背景分析参数
    do_background = args.background.lower() in ['true', 'yes', '1', 'on']

    # 获取染色体长度
    try:
        # 使用samtools view -H获取头信息，然后解析染色体长度
        cmd = ['samtools', 'view', '-H', bam_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        chrom_length = None
        for line in result.stdout.split('\n'):
            if line.startswith('@SQ') and f'SN:{target_chrom}' in line:
                for part in line.split('\t'):
                    if part.startswith('LN:'):
                        chrom_length = int(part[3:])
                        break
                break
        
        if chrom_length is None:
            print(f"[ERROR] 无法获取染色体 {target_chrom} 的长度")
            return
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] 获取染色体长度失败: {e}")
        return

    generated_files = []

    # --- 执行全长分析 ---
    if do_background:
        print(f"\n[INFO] [1/2] 正在分析 {target_chrom} 全长背景 (长度: {chrom_length/1e6:.2f} Mb)...")
        try:
            wgs_bins, wgs_counts = get_counts(bam_path, target_chrom, 0, chrom_length, wgs_bin)
            wgs_fname = os.path.join(out_dir, f"{bam_basename}_chromosome_{target_chrom}_step{wgs_bin}.png")
            plot_data(wgs_bins, wgs_counts, target_chrom, wgs_bin,
                     f"WORF-Seq Chromosome-wide Coverage\\n{target_chrom} (Step: {wgs_bin:,} bp)", 
                     wgs_fname, target_pos=target_pos)
            generated_files.append(wgs_fname)
        except Exception as e:
            print(f"[ERROR] 全染色体分析失败: {e}")
    else:
        print("[INFO] 跳过全染色体分析")

    # --- 执行精细分析 ---
    micro_bin = 500
    micro_start = max(0, target_pos - 50000)
    micro_end = min(chrom_length, target_pos + 50000)

    print(f"[INFO] [2/2] 正在分析目标区域 (+/- 50kb 范围)...")
    try:
        m_bins, m_counts = get_counts(samfile, target_chrom, micro_start, micro_end, micro_bin)
        target_fname = os.path.join(out_dir, f"{bam_basename}_target_region_{target_chrom}_{target_pos}.png")
        plot_data(m_bins, m_counts, target_chrom, micro_bin,
                 f"WORF-Seq Target Region Coverage\\n{target_chrom}:{micro_start:,}-{micro_end:,}", 
                 target_fname, target_pos=target_pos)
        generated_files.append(target_fname)
    except Exception as e:
        print(f"[ERROR] 目标区域分析失败: {e}")

    samfile.close()
    
    # 生成摘要报告
    summary_fname = os.path.join(out_dir, f"{bam_basename}_worf_seq_summary.txt")
    try:
        with open(summary_fname, 'w') as f:
            f.write("WORF-Seq Analysis Summary Report\\n")
            f.write("=" * 40 + "\\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
            f.write(f"BAM File: {bam_path}\\n")
            f.write(f"Target Chromosome: {target_chrom}\\n")
            f.write(f"Center Position: {target_pos:,}\\n")
            f.write(f"Genome-wide Step Size: {wgs_bin:,} bp\\n")
            f.write(f"Background Analysis: {'Yes' if do_background else 'No'}\\n")
            f.write(f"Chromosome Length: {chrom_length:,} bp\\n\\n")
            f.write("Generated Files:\\n")
            for file in generated_files:
                f.write(f"- {file}\\n")
        generated_files.append(summary_fname)
        print(f"[INFO] 摘要报告已保存: {summary_fname}")
    except Exception as e:
        print(f"[ERROR] 生成摘要报告失败: {e}")

    print(f"\\n[SUCCESS] 分析完成！生成了 {len(generated_files)} 个文件：")
    for file in generated_files:
        print(f"[INFO]   - {file}")

if __name__ == "__main__":
    main()