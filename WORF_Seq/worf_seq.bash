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

# 设置工作目录和日志文件
FOLDER_BASENAME=$(basename "$FOLDER_NAME")

# 检查目标目录是否有写入权限
if [[ -w "$FOLDER_NAME" ]]; then
    WORK_DIR="$FOLDER_NAME"
    LOG_FILE="${FOLDER_NAME}/${FOLDER_BASENAME}_worf_seq_pipeline.log"
    echo "[INFO] Using target directory for work: $WORK_DIR" 
else
    # 使用临时目录
    TEMP_DIR="/tmp/worf_seq_${FOLDER_BASENAME}_$$"
    mkdir -p "$TEMP_DIR"
    WORK_DIR="$TEMP_DIR"
    LOG_FILE="${TEMP_DIR}/${FOLDER_BASENAME}_worf_seq_pipeline.log"
    echo "[INFO] Target directory not writable, using temp directory: $WORK_DIR"
    echo "[WARNING] Results will be available in temp directory: $WORK_DIR"
fi

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
        echo "[ERROR] fastp not found in PATH" | tee -a "$LOG_FILE"
        exit 1
    fi
fi
echo "========================================" | tee -a "$LOG_FILE"

# 步骤2: 序列比对到参考基因组
echo "[INFO] 步骤2: 序列比对 (minimap2)" | tee -a "$LOG_FILE"
SAM_FILE="${WORK_DIR}/${FOLDER_BASENAME}_aligned_minimap.sam"
HG38_FA="${REF_DIR}/hg38.fa"
HG38_MMI="${REF_DIR}/hg38.mmi"

# 检查SAM文件是否已存在且有效
if [[ -f "$SAM_FILE" && -s "$SAM_FILE" ]]; then
    # 验证现有SAM文件的头部信息
    if head -n 1 "$SAM_FILE" | grep -q "^@"; then
        SAM_LINES=$(wc -l < "$SAM_FILE")
        echo "[SKIP] SAM文件已存在且有效，跳过序列比对步骤" | tee -a "$LOG_FILE"
        echo "[INFO] 现有文件: $SAM_FILE ($(stat -c%s "$SAM_FILE" | numfmt --to=iec)iB, $SAM_LINES lines)" | tee -a "$LOG_FILE"
    else
        echo "[WARN] 现有SAM文件无效，重新进行序列比对" | tee -a "$LOG_FILE"
        RUN_ALIGNMENT=true
    fi
else
    RUN_ALIGNMENT=true
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
        echo "[ERROR] minimap2 not found in PATH" | tee -a "$LOG_FILE"
        exit 1
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
        if python3 "$WGS_SCRIPT" \
            --bam "$BAM_FILE" \
            --chromosome "$CHROMOSOME" \
            --center "$CENTER_POSITION" \
            --step "$STEP_SIZE" \
            --background "$BACKGROUND_ANALYSIS" \
            --output "$WORK_DIR" 2>&1 | tee -a "$LOG_FILE"; then
            echo "[SUCCESS] 染色体比对图生成完成" | tee -a "$LOG_FILE"
            echo "[INFO] 生成的文件:" | tee -a "$LOG_FILE"
            for PLOT_FILE in "$EXPECTED_TARGET_PLOT" "$EXPECTED_CHROM_PLOT" "$EXPECTED_SUMMARY"; do
                if [[ -f "$PLOT_FILE" && -s "$PLOT_FILE" ]]; then
                    echo "[INFO]   - $(basename "$PLOT_FILE"): $(stat -c%s "$PLOT_FILE" | numfmt --to=iec)iB" | tee -a "$LOG_FILE"
                fi
            done
        else
            echo "[ERROR] WGSmapping.py failed" | tee -a "$LOG_FILE"
            exit 1
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
echo "[INFO]   - 输出BAM文件: $(ls -la "$BAM_FILE" 2>/dev/null | awk '{print $5}' | numfmt --to=iec)iB" | tee -a "$LOG_FILE"
echo "[INFO]   - 生成图片: $(find "$FOLDER_NAME" -name "*.png" -o -name "*.pdf" | wc -l)" | tee -a "$LOG_FILE"

echo "[SUCCESS] WORF-Seq Analysis Pipeline Completed Successfully!" | tee -a "$LOG_FILE"
echo "[INFO] Timestamp: $(date)" | tee -a "$LOG_FILE"
echo "[INFO] Log file: $LOG_FILE" | tee -a "$LOG_FILE"

# 显示结果位置信息
echo "" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "📁 RESULTS LOCATION:" | tee -a "$LOG_FILE"

if [[ "$WORK_DIR" == "$FOLDER_NAME" ]]; then
    echo "[INFO] Results are in the original directory: $WORK_DIR" | tee -a "$LOG_FILE"
    echo "[INFO] Please check the following files:" | tee -a "$LOG_FILE"
else
    echo "[WARNING] Results are in temporary directory due to permission issues:" | tee -a "$LOG_FILE"
    echo "[INFO] Temporary directory: $WORK_DIR" | tee -a "$LOG_FILE"
    echo "[INFO] Please copy results to your desired location before the system reboots" | tee -a "$LOG_FILE"
fi

echo "[INFO] Generated files:" | tee -a "$LOG_FILE"
echo "[INFO]   - BAM: $BAM_FILE" | tee -a "$LOG_FILE"
echo "[INFO]   - Target plot: $EXPECTED_TARGET_PLOT" | tee -a "$LOG_FILE"
echo "[INFO]   - Chromosome plot: $EXPECTED_CHROM_PLOT" | tee -a "$LOG_FILE"
echo "[INFO]   - Summary: $EXPECTED_SUMMARY" | tee -a "$LOG_FILE"

echo "========================================" | tee -a "$LOG_FILE"

exit 0