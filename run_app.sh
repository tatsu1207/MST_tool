#!/bin/bash
# Launch SourceTracker2 Streamlit GUI

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dada2_mst

# Check if streamlit is installed
if ! command -v streamlit &> /dev/null; then
    echo "Installing Streamlit and dependencies..."
    pip install -r "$SCRIPT_DIR/requirements.txt"
fi

# Run the app
echo "Starting SourceTracker2 GUI..."
echo "Open http://localhost:8501 in your browser"
streamlit run "$SCRIPT_DIR/app.py" --server.maxUploadSize=500
