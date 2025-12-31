# NGS Tool Analyzer

一个基于Streamlit的NGS数据分析工具，提供美观的图形化界面来运行各种生物信息学pipeline。

## 🚀 功能特性

### 🧬 已实现项目
- **Egg Indel Analysis**: CRISPR基因编辑indel分析
- **Nanobody Analysis**: 纳米抗体序列分析  

### 🚧 开发中项目
- **Evo-SEQ Analysis**: 进化序列分析 (开发中)
- **WORF-Seq Analysis**: WORF序列分析 (开发中)

## ✨ 界面特性

### 🎨 现代化设计
- **卡片式布局**: 项目选择采用响应式卡片展示
- **渐变配色**: 使用现代渐变色提升视觉效果
- **动态交互**: 悬停效果和平滑过渡动画

### 📋 智能功能
- **示例数据**: 每个项目提供完整的示例输入参数
- **文件检查**: 实时验证输入文件是否存在
- **进度监控**: 实时显示pipeline执行状态和日志
- **结果预览**: 在网页中直接查看分析结果
- **文件下载**: 支持一键下载结果文件

### 🔍 用户体验
- **错误提示**: 清晰的错误信息和解决建议
- **状态指示**: 直观的可用性状态标识
- **响应式**: 适配不同屏幕尺寸
- **日志管理**: 完整的日志展示和搜索功能
- **文件下载**: 一键下载日志和结果文件

### 📜 日志功能
- **自动保存**: Pipeline执行后自动保存日志到文件
- **实时查看**: 在网页中查看完整的执行日志
- **搜索功能**: 支持关键词搜索和高亮显示
- **统计信息**: 显示日志行数、文件大小等
- **清理工具**: 一键清理日志文件

### 📊 结果展示
- **数据预览**: CSV结果直接在网页中预览
- **统计信息**: 显示数据集基本统计指标
- **文件下载**: 支持下载结果文件和日志文件
- **压缩包内容**: 显示tar.gz文件包含的文件列表

## 📋 系统要求

- Python 3.7+
- Bash shell
- 以下生物信息学工具:
  - FLASH
  - seqkit
  - conda (可选，用于CRISPResso2环境)

## 🛠️ 安装

1. 克隆或下载项目到本地
2. 安装Python依赖:
   ```bash
   pip install -r requirements.txt
   ```
3. 确保所有生物信息学工具已安装并在PATH中

## 🎯 使用方法

### 方法1: 使用启动脚本 (推荐)

```bash
./run_streamlit.sh
```

### 方法2: 直接运行

```bash
streamlit run app.py
```

应用将在浏览器中打开，通常地址为: http://localhost:8501

## 📁 项目结构

```
NGS_Tool_syh/
├── app.py                    # Streamlit主应用
├── run_streamlit.sh          # 启动脚本
├── requirements.txt          # Python依赖
├── README.md                # 说明文档
├── Egg_Indel/               # Egg Indel分析pipeline
│   └── script/
│       └── egg_insel.bash   # Egg Indel分析脚本
├── Nanobody/                # 纳米抗体分析pipeline
│   ├── nanobody.bash        # 纳米抗体分析脚本
│   ├── trim.py              # 序列trim脚本
│   └── parse.py             # 结果解析脚本
├── Evo-SEQ/                 # Evo-SEQ分析 (开发中)
└── WORF-Seq/               # WORF-Seq分析 (开发中)
```

## 🧬 Egg Indel Analysis

**功能**: 分析CRISPR基因编辑实验中的indel突变

**输入参数**:
- 序列1文件路径 (R1)
- 序列2文件路径 (R2) 
- Barcode文件路径
- 工作名称
- 窗口大小 (可选，默认15)

**输出结果**:
- 按barcode拆分的序列文件
- CRISPResso分析结果
- 打包的结果文件 (`{工作名称}_result.tar.gz`)

## 🔬 Nanobody Analysis

**功能**: 分析纳米抗体测序数据

**输入参数**:
- 序列1文件路径 (R1)
- 序列2文件路径 (R2)
- 工作名称

**处理流程**:
1. FLASH拼接双端序列
2. 转换FASTQ为FASTA格式
3. 使用指定标记trim序列 (TGTACCTGCAGATGA...GTGACCGTGTCTTCT)
4. 解析序列并生成统计表格

**输出结果**:
- `{工作名称}_result.csv` - 分析结果表格

## 🎨 界面使用

1. **项目选择**: 在左侧边栏选择要分析的项目
2. **参数设置**: 在主界面输入所需参数
3. **执行分析**: 点击执行按钮开始分析
4. **查看结果**: 在日志区域查看执行进度和结果

## 🐛 故障排除

1. **脚本未找到错误**: 确保pipeline脚本路径正确且可执行
2. **依赖工具缺失**: 安装所需生物信息学工具 (FLASH, seqkit等)
3. **权限问题**: 确保脚本有执行权限 (`chmod +x *.bash`)

## 📝 开发说明

如需添加新的分析项目:

1. 创建对应的pipeline脚本
2. 在 `app.py` 的 `PROJECTS` 字典中添加项目配置
3. 更新 `run_script()` 函数支持新项目的参数格式

## 📄 许可证

本项目采用开源许可证，详见LICENSE文件。

## 🤝 贡献

欢迎提交Issue和Pull Request来改进这个项目！# NGS_Tool_syh
