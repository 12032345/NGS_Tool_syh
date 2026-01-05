#!/bin/bash

# WORF-Seq Analysis Pipeline
# Usage: worf_seq.bash -f folder_name -c chromosome -p center_position -s step_size -b background_analysis

# 默认参数
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="$(dirname "$SCRIPT_DIR")/WORF-Seq"
LOG_FILE=""

# 解析命令行参数
while getopts ":f:c:p:s:b:" opt; do
    case $opt in
        f) FOLDER_NAME="$OPTARG" ;;
        c) CHROMOSOME="$OPTARG" ;;
        p) CENTER_POSITION="$OPTARG" ;;
        s) STEP_SIZE="$OPTARG" ;;
        b) BACKGROUND_ANALYSIS="$OPTARG" ;;
        \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done

# 检查必需参数
if [[ -z "$FOLDER_NAME" || -z "$CHROMOSOME" || -z "$CENTER_POSITION" ]]; then
    echo "[ERROR] Missing required parameters"
    echo "Usage: $0 -f folder_name -c chromosome -p center_position [-s step_size] [-b background_analysis]"
    exit 1
fi

# 设置默认值
STEP_SIZE=${STEP_SIZE:-100000}
BACKGROUND_ANALYSIS=${BACKGROUND_ANALYSIS:-true}

# 设置工作目录和日志文件：统一使用时间戳临时目录，避免目标目录权限问题
FOLDER_BASENAME=$(basename "$FOLDER_NAME")
TIMESTAMP=$(date +%s)
WORK_DIR="/tmp/worf_seq_${FOLDER_BASENAME}_${TIMESTAMP}"
mkdir -p "$WORK_DIR"
LOG_FILE="${WORK_DIR}/${FOLDER_BASENAME}_worf_seq_pipeline.log"
echo "[INFO] Using temp work directory: $WORK_DIR" 
echo "[INFO] 原始数据仍从目录读取: $FOLDER_NAME"

# 确保管道失败能正确上报（防止 tee 掩盖子进程退出码）
set -o pipefail

echo "[INFO] WORF-Seq Analysis Pipeline Started" | tee "$LOG_FILE"
echo "[INFO] Timestamp: $(date)" | tee -a "$LOG_FILE"
echo "[INFO] Parameters:" | tee -a "$LOG_FILE"
echo "[INFO]   - Folder: $FOLDER_NAME" | tee -a "$LOG_FILE"
echo "[INFO]   - Chromosome: $CHROMOSOME" | tee -a "$LOG_FILE"
echo "[INFO]   - Center Position: $CENTER_POSITION" | tee -a "$LOG_FILE"
echo "[INFO]   - Step Size: $STEP_SIZE" | tee -a "$LOG_FILE"
echo "[INFO]   - Background Analysis: $BACKGROUND_ANALYSIS" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

# 检查输入文件是否存在
FOLDER_BASENAME=$(basename "$FOLDER_NAME")
RAW_R1="${FOLDER_NAME}/${FOLDER_BASENAME}_raw_1.fq.gz"
RAW_R2="${FOLDER_NAME}/${FOLDER_BASENAME}_raw_2.fq.gz"

if [[ ! -f "$RAW_R1" ]]; then
    echo "[ERROR] Raw file not found: $RAW_R1" | tee -a "$LOG_FILE"
    exit 1
fi

if [[ ! -f "$RAW_R2" ]]; then
    echo "[ERROR] Raw file not found: $RAW_R2" | tee -a "$LOG_FILE"
    exit 1
fi

echo "[INFO] Input files verified:" | tee -a "$LOG_FILE"
echo "[INFO]   - R1: $RAW_R1" | tee -a "$LOG_FILE"
echo "[INFO]   - R2: $RAW_R2" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

# 步骤1: 质控处理
echo "[INFO] 步骤1: 开始质控处理 (fastp)" | tee -a "$LOG_FILE"
CLEAN_R1="${WORK_DIR}/${FOLDER_BASENAME}_clean_1.fq.gz"
CLEAN_R2="${WORK_DIR}/${FOLDER_BASENAME}_clean_2.fq.gz"

# 检查质控结果文件是否已存在
if [[ -f "$CLEAN_R1" && -f "$CLEAN_R2" && -s "$CLEAN_R1" && -s "$CLEAN_R2" ]]; then
    echo "[SKIP] 质控文件已存在，跳过质控步骤" | tee -a "$LOG_FILE"
    echo "[INFO] 现有文件:" | tee -a "$LOG_FILE"
    echo "[INFO]   - Clean R1: $CLEAN_R1 ($(stat -c%s "$CLEAN_R1" | numfmt --to=iec)iB)" | tee -a "$LOG_FILE"
    echo "[INFO]   - Clean R2: $CLEAN_R2 ($(stat -c%s "$CLEAN_R2" | numfmt --to=iec)iB)" | tee -a "$LOG_FILE"
else
    if command -v fastp >/dev/null 2>&1; then
        echo "[INFO] Running fastp..." | tee -a "$LOG_FILE"
        if fastp -i "$RAW_R1" -I "$RAW_R2" -o "$CLEAN_R1" -O "$CLEAN_R2" -Q -L 2>&1 | tee -a "$LOG_FILE"; then
            echo "[SUCCESS] 质控完成" | tee -a "$LOG_FILE"
            echo "[INFO] Clean files:" | tee -a "$LOG_FILE"
            echo "[INFO]   - Clean R1: $CLEAN_R1 ($(stat -c%s "$CLEAN_R1" | numfmt --to=iec)iB)" | tee -a "$LOG_FILE"
            echo "[INFO]   - Clean R2: $CLEAN_R2 ($(stat -c%s "$CLEAN_R2" | numfmt --to=iec)iB)" | tee -a "$LOG_FILE"
        else
            echo "[ERROR] fastp failed" | tee -a "$LOG_FILE"
            exit 1
        fi
    else
        # fastp 缺失时兜底：若已存在清洗文件则继续，否则使用原始文件作为清洗输入
        if [[ -f "$CLEAN_R1" && -f "$CLEAN_R2" && -s "$CLEAN_R1" && -s "$CLEAN_R2" ]]; then
            echo "[WARN] fastp 未安装，但检测到现有clean文件，继续后续流程" | tee -a "$LOG_FILE"
        else
            echo "[WARN] fastp 未安装，使用原始文件直接作为 clean 输入" | tee -a "$LOG_FILE"
            ln -sf "$RAW_R1" "$CLEAN_R1"
            ln -sf "$RAW_R2" "$CLEAN_R2"
            echo "[INFO] 已创建符号链接替代 clean 文件" | tee -a "$LOG_FILE"
        fi
    fi
fi
echo "========================================" | tee -a "$LOG_FILE"

# 步骤2: 序列比对到参考基因组
echo "[INFO] 步骤2: 序列比对 (minimap2)" | tee -a "$LOG_FILE"
SAM_FILE="${WORK_DIR}/${FOLDER_BASENAME}_aligned_minimap.sam"
# 提前定义BAM路径用于跳过比对的检测
BAM_FILE="${WORK_DIR}/${FOLDER_BASENAME}_aligned_minimap.sorted.bam"
BAM_INDEX="${BAM_FILE}.bai"
HG38_FA="${REF_DIR}/hg38.fa"
HG38_MMI="${REF_DIR}/hg38.mmi"
RUN_ALIGNMENT=true

# 如果之前因权限使用过临时目录，尝试复用已有 SAM/BAM（避免重复跑 minimap2）
TMP_GLOB="/tmp/worf_seq_${FOLDER_BASENAME}_*"
FOUND_TMP_BAM=$(ls $TMP_GLOB/${FOLDER_BASENAME}_aligned_minimap.sorted.bam 2>/dev/null | head -n 1)
FOUND_TMP_SAM=$(ls $TMP_GLOB/${FOLDER_BASENAME}_aligned_minimap.sam 2>/dev/null | head -n 1)

if [[ ! -s "$BAM_FILE" && -n "$FOUND_TMP_BAM" ]]; then
    echo "[INFO] 复用临时目录中的BAM: $FOUND_TMP_BAM" | tee -a "$LOG_FILE"
    ln -sf "$FOUND_TMP_BAM" "$BAM_FILE"
    if [[ -f "${FOUND_TMP_BAM}.bai" ]]; then
        ln -sf "${FOUND_TMP_BAM}.bai" "$BAM_INDEX"
    fi
fi

if [[ ! -s "$SAM_FILE" && -n "$FOUND_TMP_SAM" ]]; then
    echo "[INFO] 复用临时目录中的SAM: $FOUND_TMP_SAM" | tee -a "$LOG_FILE"
    ln -sf "$FOUND_TMP_SAM" "$SAM_FILE"
fi

# 若已存在排序BAM，直接跳过比对
if [[ -f "$BAM_FILE" && -s "$BAM_FILE" ]]; then
    echo "[SKIP] 检测到已存在排序BAM，跳过序列比对" | tee -a "$LOG_FILE"
    RUN_ALIGNMENT=false
else
    # 检查SAM文件是否已存在且有效
    if [[ -f "$SAM_FILE" && -s "$SAM_FILE" ]]; then
        # 验证现有SAM文件的头部信息
        if head -n 1 "$SAM_FILE" | grep -q "^@"; then
            SAM_LINES=$(wc -l < "$SAM_FILE")
            echo "[SKIP] SAM文件已存在且有效，跳过序列比对步骤" | tee -a "$LOG_FILE"
            echo "[INFO] 现有文件: $SAM_FILE ($(stat -c%s "$SAM_FILE" | numfmt --to=iec)iB, $SAM_LINES lines)" | tee -a "$LOG_FILE"
            RUN_ALIGNMENT=false
        else
            echo "[WARN] 现有SAM文件无效，重新进行序列比对" | tee -a "$LOG_FILE"
            RUN_ALIGNMENT=true
        fi
    else
        RUN_ALIGNMENT=true
    fi
fi

if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    if command -v minimap2 >/dev/null 2>&1; then
        echo "[INFO] Running minimap2..." | tee -a "$LOG_FILE"
        # 优先使用预构建的索引文件，如果不存在则使用FASTA文件
        if [[ -f "$HG38_MMI" ]]; then
            echo "[INFO] Using pre-built index: $HG38_MMI" | tee -a "$LOG_FILE"
            if minimap2 -ax sr -t 8 "$HG38_MMI" "$CLEAN_R1" "$CLEAN_R2" > "$SAM_FILE" 2>> "$LOG_FILE"; then
                echo "[SUCCESS] 序列比对完成" | tee -a "$LOG_FILE"
                echo "[INFO] SAM file: $SAM_FILE" | tee -a "$LOG_FILE"
            else
                echo "[ERROR] minimap2 failed" | tee -a "$LOG_FILE"
                exit 1
            fi
        elif [[ -f "$HG38_FA" ]]; then
            echo "[INFO] Using reference genome: $HG38_FA" | tee -a "$LOG_FILE"
            if minimap2 -ax sr -t 8 "$HG38_FA" "$CLEAN_R1" "$CLEAN_R2" > "$SAM_FILE" 2>> "$LOG_FILE"; then
                echo "[SUCCESS] 序列比对完成" | tee -a "$LOG_FILE"
                echo "[INFO] SAM file: $SAM_FILE" | tee -a "$LOG_FILE"
            else
                echo "[ERROR] minimap2 failed" | tee -a "$LOG_FILE"
                exit 1
            fi
        else
            echo "[ERROR] Reference genome files not found:" | tee -a "$LOG_FILE"
            echo "[ERROR]   - Index file: $HG38_MMI" | tee -a "$LOG_FILE"
            echo "[ERROR]   - FASTA file: $HG38_FA" | tee -a "$LOG_FILE"
            exit 1
        fi
    else
        # minimap2 缺失时兜底：若已有有效 SAM/BAM 则跳过；否则退出
        if [[ -f "$BAM_FILE" && -s "$BAM_FILE" ]]; then
            echo "[WARN] minimap2 未安装，检测到现有BAM，跳过比对步骤" | tee -a "$LOG_FILE"
            RUN_ALIGNMENT=false
        elif [[ -f "$SAM_FILE" && -s "$SAM_FILE" ]]; then
            echo "[WARN] minimap2 未安装，检测到现有SAM，跳过比对步骤" | tee -a "$LOG_FILE"
            RUN_ALIGNMENT=false
        else
            echo "[ERROR] minimap2 not found in PATH，且未检测到已有SAM/BAM可用" | tee -a "$LOG_FILE"
            exit 1
        fi
    fi
fi
echo "========================================" | tee -a "$LOG_FILE"


echo "========================================" | tee -a "$LOG_FILE"

# 步骤3: SAM转换为BAM文件
echo "[INFO] 步骤3: SAM转BAM (samtools)" | tee -a "$LOG_FILE"
BAM_FILE="${WORK_DIR}/${FOLDER_BASENAME}_aligned_minimap.sorted.bam"
BAM_INDEX="${BAM_FILE}.bai"

# 检查BAM文件和索引是否已存在
if [[ -f "$BAM_FILE" && -s "$BAM_FILE" && -f "$BAM_INDEX" && -s "$BAM_INDEX" ]]; then
    echo "[SKIP] BAM文件和索引已存在，跳过SAM转BAM步骤" | tee -a "$LOG_FILE"
    echo "[INFO] 现有文件:" | tee -a "$LOG_FILE"
    echo "[INFO]   - BAM file: $BAM_FILE ($(stat -c%s "$BAM_FILE" | numfmt --to=iec)iB)" | tee -a "$LOG_FILE"
    echo "[INFO]   - BAM index: $BAM_INDEX ($(stat -c%s "$BAM_INDEX" | numfmt --to=iec)iB)" | tee -a "$LOG_FILE"
else
    if command -v samtools >/dev/null 2>&1; then
        # 检查是否需要排序
        if [[ -f "$BAM_FILE" && -s "$BAM_FILE" ]]; then
            echo "[INFO] BAM文件已存在，检查索引..." | tee -a "$LOG_FILE"
            NEED_INDEX=true
        else
            echo "[INFO] Sorting SAM file..." | tee -a "$LOG_FILE"
            if samtools sort -@ 8 -o "$BAM_FILE" "$SAM_FILE" 2>&1 | tee -a "$LOG_FILE"; then
                echo "[SUCCESS] BAM sorting completed" | tee -a "$LOG_FILE"
                echo "[INFO] BAM file: $BAM_FILE ($(stat -c%s "$BAM_FILE" | numfmt --to=iec)iB)" | tee -a "$LOG_FILE"
                NEED_INDEX=true
            else
                echo "[ERROR] BAM sorting failed" | tee -a "$LOG_FILE"
                exit 1
            fi
        fi
        
        if [[ "$NEED_INDEX" == "true" ]]; then
            echo "[INFO] Indexing BAM file..." | tee -a "$LOG_FILE"
            if samtools index "$BAM_FILE" 2>&1 | tee -a "$LOG_FILE"; then
                echo "[SUCCESS] BAM indexing completed" | tee -a "$LOG_FILE"
                echo "[INFO] BAM index: $BAM_INDEX ($(stat -c%s "$BAM_INDEX" | numfmt --to=iec)iB)" | tee -a "$LOG_FILE"
            else
                echo "[ERROR] BAM indexing failed" | tee -a "$LOG_FILE"
                exit 1
            fi
        fi
    else
        echo "[ERROR] samtools not found in PATH" | tee -a "$LOG_FILE"
        exit 1
    fi
fi
echo "========================================" | tee -a "$LOG_FILE"

# 步骤4: 染色体比对图生成
echo "[INFO] 步骤4: 染色体比对图生成 (WGSmapping.py)" | tee -a "$LOG_FILE"
WGS_SCRIPT="${SCRIPT_DIR}/WGSmapping.py"

# 预期的输出文件名
FOLDER_BASENAME=$(basename "$FOLDER_NAME")
EXPECTED_TARGET_PLOT="${WORK_DIR}/${FOLDER_BASENAME}_target_region_${CHROMOSOME}_${CENTER_POSITION}.png"
EXPECTED_CHROM_PLOT="${WORK_DIR}/${FOLDER_BASENAME}_chromosome_${CHROMOSOME}_step${STEP_SIZE}.png"
EXPECTED_SUMMARY="${WORK_DIR}/${FOLDER_BASENAME}_worf_seq_summary.txt"

# 检查图表文件是否已存在
PLOTS_EXIST=true
for PLOT_FILE in "$EXPECTED_TARGET_PLOT" "$EXPECTED_CHROM_PLOT" "$EXPECTED_SUMMARY"; do
    if [[ "$BACKGROUND_ANALYSIS" == "true" ]]; then
        # 如果背景分析为true，检查所有文件
        if [[ ! -f "$PLOT_FILE" || ! -s "$PLOT_FILE" ]]; then
            PLOTS_EXIST=false
            break
        fi
    else
        # 如果背景分析为false，只检查目标区域图和摘要
        if [[ "$PLOT_FILE" == "$EXPECTED_CHROM_PLOT" ]]; then
            continue
        fi
        if [[ ! -f "$PLOT_FILE" || ! -s "$PLOT_FILE" ]]; then
            PLOTS_EXIST=false
            break
        fi
    fi
done

if [[ "$PLOTS_EXIST" == "true" ]]; then
    echo "[SKIP] 图表文件已存在，跳过染色体比对图生成步骤" | tee -a "$LOG_FILE"
    echo "[INFO] 现有文件:" | tee -a "$LOG_FILE"
    for PLOT_FILE in "$EXPECTED_TARGET_PLOT" "$EXPECTED_CHROM_PLOT" "$EXPECTED_SUMMARY"; do
        if [[ -f "$PLOT_FILE" && -s "$PLOT_FILE" ]]; then
            echo "[INFO]   - $(basename "$PLOT_FILE"): $(stat -c%s "$PLOT_FILE" | numfmt --to=iec)iB" | tee -a "$LOG_FILE"
        fi
    done
else
    if [[ -f "$WGS_SCRIPT" ]]; then
        echo "[INFO] Running WGSmapping.py..." | tee -a "$LOG_FILE"
        # 运行 WGSmapping 并捕获子进程的退出码（避免 tee 掩盖）
        python3 "$WGS_SCRIPT" \
            --bam "$BAM_FILE" \
            --chromosome "$CHROMOSOME" \
            --center "$CENTER_POSITION" \
            --step "$STEP_SIZE" \
            --background "$BACKGROUND_ANALYSIS" \
            --output "$WORK_DIR" 2>&1 | tee -a "$LOG_FILE"
        PY_EXIT=${PIPESTATUS[0]}
        if [[ $PY_EXIT -eq 0 ]]; then
            echo "[SUCCESS] 染色体比对图生成完成" | tee -a "$LOG_FILE"
            echo "[INFO] 生成的文件:" | tee -a "$LOG_FILE"
            for PLOT_FILE in "$EXPECTED_TARGET_PLOT" "$EXPECTED_CHROM_PLOT" "$EXPECTED_SUMMARY"; do
                if [[ -f "$PLOT_FILE" && -s "$PLOT_FILE" ]]; then
                    echo "[INFO]   - $(basename "$PLOT_FILE"): $(stat -c%s "$PLOT_FILE" | numfmt --to=iec)iB" | tee -a "$LOG_FILE"
                fi
            done
        else
            echo "[ERROR] WGSmapping.py failed (exit code $PY_EXIT)" | tee -a "$LOG_FILE"
            exit $PY_EXIT
        fi
    else
        echo "[ERROR] WGSmapping.py not found: $WGS_SCRIPT" | tee -a "$LOG_FILE"
        exit 1
    fi
fi
echo "========================================" | tee -a "$LOG_FILE"

# 统计信息
echo "[INFO] 分析完成统计:" | tee -a "$LOG_FILE"
echo "[INFO]   - 输入文件: 2" | tee -a "$LOG_FILE"
if [[ -f "$BAM_FILE" && -s "$BAM_FILE" ]]; then
    echo "[INFO]   - 输出BAM文件: $(stat -c%s "$BAM_FILE" | numfmt --to=iec)iB" | tee -a "$LOG_FILE"
else
    echo "[INFO]   - 输出BAM文件: (not found)" | tee -a "$LOG_FILE"
fi
IMG_COUNT=$(find "$WORK_DIR" -maxdepth 1 -type f \( -name "*.png" -o -name "*.pdf" \) | wc -l)
echo "[INFO]   - 生成图片: $IMG_COUNT" | tee -a "$LOG_FILE"

echo "[SUCCESS] WORF-Seq Analysis Pipeline Completed Successfully!" | tee -a "$LOG_FILE"
echo "[INFO] Timestamp: $(date)" | tee -a "$LOG_FILE"
echo "[INFO] Log file: $LOG_FILE" | tee -a "$LOG_FILE"

# 显示结果位置信息
echo "" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "📁 RESULTS LOCATION:" | tee -a "$LOG_FILE"
echo "[INFO] Results are in temporary directory: $WORK_DIR" | tee -a "$LOG_FILE"
echo "[INFO] Please copy results to your desired location before the system reboots" | tee -a "$LOG_FILE"

echo "[INFO] Generated files:" | tee -a "$LOG_FILE"
if [[ -f "$BAM_FILE" && -s "$BAM_FILE" ]]; then
    echo "[INFO]   - BAM: $BAM_FILE" | tee -a "$LOG_FILE"
fi
if [[ -f "$EXPECTED_TARGET_PLOT" && -s "$EXPECTED_TARGET_PLOT" ]]; then
    echo "[INFO]   - Target plot: $EXPECTED_TARGET_PLOT" | tee -a "$LOG_FILE"
fi
if [[ -f "$EXPECTED_CHROM_PLOT" && -s "$EXPECTED_CHROM_PLOT" ]]; then
    echo "[INFO]   - Chromosome plot: $EXPECTED_CHROM_PLOT" | tee -a "$LOG_FILE"
fi
if [[ -f "$EXPECTED_SUMMARY" && -s "$EXPECTED_SUMMARY" ]]; then
    echo "[INFO]   - Summary: $EXPECTED_SUMMARY" | tee -a "$LOG_FILE"
fi

echo "========================================" | tee -a "$LOG_FILE"

exit 0