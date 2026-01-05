#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="ngs_tools"
ENV_FILE="$PROJECT_DIR/environment.yml"

echo "🚀 NGS Tool Analyzer 启动脚本"
echo "================================"

# 如果 conda 可用，优先使用 conda run 来保证在指定环境中执行
if command -v conda >/dev/null 2>&1; then
    # 检查环境是否存在
    if conda env list | awk '{print $1}' | grep -q "^${ENV_NAME}$"; then
        echo "🔁 已检测到 conda 环境: ${ENV_NAME}，将使用该环境运行 Streamlit"
    else
        # 如果 environment.yml 存在，则创建环境
        if [ -f "$ENV_FILE" ]; then
            echo "🛠 未检测到 ${ENV_NAME} 环境，正在根据 environment.yml 创建..."
            conda env create -f "$ENV_FILE" -n "$ENV_NAME" || {
                echo "❌ 创建环境失败，请手动检查 environment.yml"; exit 1;
            }
            echo "✅ 环境 ${ENV_NAME} 已创建"
        else
            echo "❌ 未找到 environment.yml，无法自动创建 ${ENV_NAME} 环境"; exit 1
        fi
    fi

    # 确认 app.py 存在
    if [ ! -f "$PROJECT_DIR/app.py" ]; then
        echo "❌ 未找到 app.py，请在项目根目录运行此脚本"; exit 1
    fi

    echo "📦 使用 conda 环境: $ENV_NAME 运行 Streamlit"
    echo "🌐 应用将在浏览器中打开 (默认端口 8501 或第一个空闲端口)"
    echo "⏹️ 按 Ctrl+C 停止应用"

    # 查找一个空闲端口（优先 8501）
    PORT=8501
    while ss -ltn | awk '{print $4}' | grep -q ":${PORT}\$"; do
        PORT=$((PORT+1))
    done

    echo "➡️ 使用端口: $PORT"

    # 使用 conda run 启动 Streamlit，保持在前台
    exec conda run -n "$ENV_NAME" streamlit run "$PROJECT_DIR/app.py" --server.port "$PORT" --server.address 0.0.0.0 --server.headless true
else
    echo "⚠️ conda 未安装，尝试在当前 Python 环境中直接运行 Streamlit"
    if ! command -v streamlit >/dev/null 2>&1; then
        echo "❌ 未找到 streamlit，请先安装或安装 conda 并创建 ${ENV_NAME} 环境"
        exit 1
    fi
    exec streamlit run "$PROJECT_DIR/app.py"
fi