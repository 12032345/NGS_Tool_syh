#!/usr/bin/env python3
"""
内置 barcode 数据管理
包含96个标准 barcode 序列
"""

import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def load_barcodes_from_file(file_path):
    """从文件中加载 barcode 序列"""
    barcodes = {}
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # 过滤空行并去除空白字符
        sequences = [line.strip() for line in lines if line.strip()]
        
        # 为每个序列分配序号（1-96）
        for i, sequence in enumerate(sequences, 1):
            if i <= 96:
                barcodes[i] = sequence
            else:
                break
                
        return barcodes
    except Exception as e:
        print(f"加载 barcode 文件失败: {e}")
        # 如果加载失败，返回默认的 barcode
        return get_default_barcodes()

def get_default_barcodes():
    """获取默认的 barcode 序列（备用）"""
    return {
        1: "ATCACGTT",
        2: "CGATGTAT", 
        3: "TTAGGCAT",
        4: "TGACCTCA",
        5: "ACAGTGGA",
        6: "GCCAATGT"
    }

# 从文件加载 barcode 序列（项目内相对路径）
BARCODE_FILE = os.path.join(BASE_DIR, "Egg_Indel", "barcode.txt")
BARCODES = load_barcodes_from_file(BARCODE_FILE)

def get_barcode_sequence(barcode_num):
    """根据 barcode 序号获取序列"""
    return BARCODES.get(barcode_num, "")

def generate_barcode_file(barcode_num, output_path):
    """生成包含指定 barcode 的文件"""
    sequence = get_barcode_sequence(barcode_num)
    if not sequence:
        return False, f"无效的 barcode 序号: {barcode_num}"
    
    try:
        with open(output_path, 'w') as f:
            f.write(f"barcode_{barcode_num:02d}\t{sequence}\n")
        return True, f"成功生成 {barcode_num} 号 barcode 文件"
    except Exception as e:
        return False, f"生成 barcode 文件失败: {e}"

def list_all_barcodes():
    """列出所有 barcode"""
    return [(i, seq) for i, seq in BARCODES.items()]

def get_barcode_display_name(barcode_num):
    """获取 barcode 的显示名称"""
    return f"Barcode {barcode_num:02d} ({BARCODES.get(barcode_num, 'Unknown')})"