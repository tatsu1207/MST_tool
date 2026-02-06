#!/bin/bash
# Launch SourceTracker2 Streamlit GUI

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Find and initialize conda
init_conda() {
    # Common conda installation paths
    local conda_paths=(
        "$HOME/miniforge"
        "$HOME/miniforge3"
        "$HOME/mambaforge"
        "$HOME/miniconda3"
        "$HOME/miniconda"
        "$HOME/anaconda3"
        "/opt/conda"
        "/opt/miniforge3"
    )

    for conda_base in "${conda_paths[@]}"; do
        if [ -f "$conda_base/etc/profile.d/conda.sh" ]; then
            source "$conda_base/etc/profile.d/conda.sh"
            return 0
        fi
    done

    # Try using conda info if conda is in PATH
    if command -v conda &> /dev/null; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        return 0
    fi

    echo "Error: Could not find conda installation."
    echo "Please ensure Miniforge/Miniconda is installed."
    return 1
}

# Initialize conda
if ! init_conda; then
    exit 1
fi

# Activate conda environment
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
