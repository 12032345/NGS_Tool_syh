import pysam
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime
import argparse

def get_counts(samfile, chrom, start, end, bin_size):
    """计算指定区间内每个bin的符合条件的reads数"""
    bins = range(start, end, bin_size)
    counts = []
    
    for b_start in bins:
        b_end = b_start + bin_size
        bin_count = 0
        # fetch 利用索引快速定位
        try:
            for read in samfile.fetch(chrom, b_start, b_end):
                if read.is_unmapped or read.is_secondary:
                    continue
                
                # 计算 read 在当前 bin 内的覆盖长度
                overlap = read.get_overlap(b_start, b_end)
                read_len = read.query_length if read.query_length > 0 else 1
                
                # 核心算法：重叠部分 > 50%
                if (overlap / read_len) > 0.5:
                    bin_count += 1
        except ValueError:
            # 如果染色体名称在BAM里找不到，跳过
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
    parser.add_argument('--bam', help='排序后的 BAM 文件路径', required=False)
    parser.add_argument('--chrom', help='目标染色体 (例如 chr6)', required=False)
    parser.add_argument('--pos', help='目标中心位置 (整数)', type=int, required=False)
    parser.add_argument('--wgs-bin', help='全染色体绘图步长 (bp)', type=int, default=100000)
    parser.add_argument('--skip-wgs', help='跳过全染色体分析，只做目标微区分析', action='store_true')

    args = parser.parse_args()

    print("=== 捕获测序富集分析工具 ===")

    # 1. 获取 BAM 路径
    if args.bam:
        bam_path = args.bam
    else:
        while True:
            bam_path = input("请输入服务器中排序后的 BAM 文件路径: ").strip()
            if os.path.exists(bam_path):
                if os.path.exists(bam_path + ".bai"):
                    break
                else:
                    print("❌ 警告：未找到索引文件 (.bai)。请先运行 'samtools index'。")
            else:
                print("❌ 错误：文件路径不存在，请重新输入。")

    try:
        samfile = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"❌ 无法读取BAM文件: {e}")
        return

    # 输出目录：在 BAM 所在目录下创建 WGSmapping 子目录
    bam_abspath = os.path.abspath(bam_path)
    bam_dir = os.path.dirname(bam_abspath) or os.getcwd()
    out_dir = os.path.join(bam_dir, 'WGSmapping')
    os.makedirs(out_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # 2. 全长绘图间隔
    wgs_bin = args.wgs_bin

    # 3. 关注的基因位置
    if args.chrom and args.pos:
        target_chrom = args.chrom
        target_pos = args.pos
    else:
        print("\n2. 请输入目标区域信息：")
        target_chrom = input("   染色体编号 (如 chr3): ").strip()
        target_pos = int(input("   中心点位置 (如 150176116): ").strip())

    chrom_length = samfile.get_reference_length(target_chrom)

    # --- 执行全长分析 ---
    if not args.skip_wgs:
        print(f"\n[1/2] 正在分析 {target_chrom} 全长背景 (长度: {chrom_length/1e6:.2f} Mb)...")
        wgs_bins, wgs_counts = get_counts(samfile, target_chrom, 0, chrom_length, wgs_bin)
        wgs_fname = os.path.join(out_dir, f"WGS_Overview_{target_chrom}_{timestamp}.png")
        plot_data(wgs_bins, wgs_counts, target_chrom, wgs_bin,
              f"WGS Background: {target_chrom}", wgs_fname, target_pos=target_pos)
    else:
        print("跳过全染色体分析（--skip-wgs 指定）")

    # --- 执行精细分析 ---
    micro_bin = 500
    micro_start = max(0, target_pos - 50000)
    micro_end = min(chrom_length, target_pos + 50000)

    print(f"[2/2] 正在分析目标区域 (+/- 50kb 范围)...")
    m_bins, m_counts = get_counts(samfile, target_chrom, micro_start, micro_end, micro_bin)
    target_fname = os.path.join(out_dir, f"Target_Detail_{target_chrom}_{timestamp}.png")
    plot_data(m_bins, m_counts, target_chrom, micro_bin,
              f"Target Enrichment: {target_chrom}:{target_pos}", target_fname, target_pos=target_pos)

    samfile.close()
    print("\n🎉 分析完成！请检查当前目录下的 PNG 图片文件。")

if __name__ == "__main__":
    main()
