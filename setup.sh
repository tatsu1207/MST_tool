#!/bin/bash

# DADA2 + SourceTracker2 Environment Setup
# Creates conda environment "dada2_mst" with all required dependencies.

# Color definitions
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

LOG_FILE="dada2_mst_setup.log"
ENV_NAME="dada2_mst"

# Clear previous log
echo "--- DADA2 MST Setup Log: $(date) ---" > "$LOG_FILE"

echo -e "${BLUE}===============================================${NC}"
echo -e "${BLUE}   DADA2 + SourceTracker2 Environment Setup   ${NC}"
echo -e "${BLUE}   Environment name: ${ENV_NAME}              ${NC}"
echo -e "${BLUE}===============================================${NC}"
echo -e "${YELLOW}Logs: $LOG_FILE${NC}"

# Function to check command success
check_status() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ $1${NC}"
    else
        echo -e "${RED}✗ $1 failed. Check $LOG_FILE${NC}"
        exit 1
    fi
}

# 1. Check for Conda/Mamba
echo -e "\n${BLUE}[1/5] Checking Dependencies...${NC}"

if ! command -v conda &> /dev/null; then
    # Try to source conda from common install locations
    for CONDA_SH in "$HOME/miniforge3/etc/profile.d/conda.sh" \
                     "$HOME/miniconda3/etc/profile.d/conda.sh" \
                     "$HOME/anaconda3/etc/profile.d/conda.sh"; do
        if [ -f "$CONDA_SH" ]; then
            source "$CONDA_SH"
            break
        fi
    done
fi

if ! command -v conda &> /dev/null; then
    echo -e "${RED}✗ Conda not found. Install Miniforge first:${NC}"
    echo -e "${YELLOW}wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh${NC}"
    echo -e "${YELLOW}bash Miniforge3-Linux-x86_64.sh${NC}"
    exit 1
fi

# Initialize conda
eval "$(conda shell.bash hook)"

# Install Mamba if needed
if ! command -v mamba &> /dev/null; then
    echo -e "${YELLOW}Installing Mamba...${NC}"
    conda install -n base -c conda-forge mamba -y &>> "$LOG_FILE"
    check_status "Mamba installation"
else
    echo -e "${GREEN}✓ Mamba already installed${NC}"
fi

# 2. Install system dependencies (optional, for WSL2/Linux)
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo -e "\n${BLUE}[2/5] Installing system libraries...${NC}"
    echo -e "${YELLOW}(Requires sudo - you may be prompted for password)${NC}"

    sudo apt-get update &>> "$LOG_FILE"
    sudo apt-get install -y git build-essential libglpk-dev libxml2-dev \
        libssl-dev libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev \
        libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        gfortran libgfortran5 &>> "$LOG_FILE"
    check_status "System libraries"
else
    echo -e "\n${BLUE}[2/5] Skipping system libraries (not Linux)${NC}"
fi

# 3. Configure channels
echo -e "\n${BLUE}[3/5] Configuring conda channels...${NC}"
conda config --set always_yes yes &>> "$LOG_FILE"
conda config --add channels defaults &>> "$LOG_FILE"
conda config --add channels bioconda &>> "$LOG_FILE"
conda config --add channels conda-forge &>> "$LOG_FILE"
conda config --set channel_priority flexible &>> "$LOG_FILE"
check_status "Channel configuration"

# 4. Remove existing environment if present
if conda env list | grep -q "^$ENV_NAME "; then
    echo -e "\n${YELLOW}Environment '$ENV_NAME' exists.${NC}"
    read -p "Remove and recreate? (y/n): " confirm
    if [[ $confirm == [yY] ]]; then
        echo -e "${BLUE}Removing environment...${NC}"
        mamba env remove -n $ENV_NAME -y &>> "$LOG_FILE"
        check_status "Environment removal"
    else
        echo -e "${YELLOW}Keeping existing environment. Exiting.${NC}"
        exit 0
    fi
fi

# 5. Create environment
echo -e "\n${BLUE}[4/5] Creating environment...${NC}"
echo -e "${YELLOW}This may take 10-15 minutes...${NC}"

mamba create -n $ENV_NAME \
    python=3.9 \
    r-base=4.2 \
    r-essentials \
    bioconductor-dada2 \
    r-optparse \
    cutadapt \
    seqkit \
    biom-format \
    h5py \
    scikit-bio \
    numpy \
    scipy \
    pandas \
    cython \
    icu=58 \
    libstdcxx-ng \
    -c conda-forge -c bioconda -c defaults -y &>> "$LOG_FILE"

check_status "Environment creation"

# Activate environment
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

# Set up LD_LIBRARY_PATH (for Linux/WSL2)
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    mkdir -p $CONDA_PREFIX/etc/conda/activate.d
    mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d

    cat > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh << 'EOF'
#!/bin/bash
export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
EOF

    cat > $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh << 'EOF'
#!/bin/bash
export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
unset OLD_LD_LIBRARY_PATH
EOF

    # Re-activate to apply env vars
    conda deactivate
    conda activate $ENV_NAME
fi

# 6. Install pip packages
echo -e "\n${BLUE}[5/5] Installing SourceTracker2 and visualization packages...${NC}"

pip install --upgrade pip &>> "$LOG_FILE"
check_status "Pip upgrade"

echo -e "${YELLOW}Installing SourceTracker2...${NC}"
pip install git+https://github.com/biom-format/biom-format.git &>> "$LOG_FILE"
pip install git+https://github.com/biocore/sourcetracker2.git &>> "$LOG_FILE"
check_status "SourceTracker2"

echo -e "${YELLOW}Installing visualization packages...${NC}"
pip install streamlit plotly matplotlib seaborn biopython &>> "$LOG_FILE"
check_status "Visualization packages"

# 7. Verification
echo -e "\n${BLUE}===============================================${NC}"
echo -e "${BLUE}           VERIFYING INSTALLATION              ${NC}"
echo -e "${BLUE}===============================================${NC}"

# Test SourceTracker2
python -c "import sourcetracker; print('SourceTracker2:', sourcetracker.__version__)" &>> "$LOG_FILE"
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ SourceTracker2: INSTALLED${NC}"
else
    echo -e "${YELLOW}! SourceTracker2: Check manually${NC}"
fi

# Test Cutadapt
if command -v cutadapt &> /dev/null; then
    VERSION=$(cutadapt --version 2>&1)
    echo -e "${GREEN}✓ Cutadapt: $VERSION${NC}"
else
    echo -e "${YELLOW}! Cutadapt: Not found${NC}"
fi

# Test DADA2
Rscript -e "library(dada2); cat('DADA2:', as.character(packageVersion('dada2')), '\n')" &>> "$LOG_FILE"
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ DADA2: INSTALLED${NC}"
else
    echo -e "${YELLOW}! DADA2: Check manually${NC}"
fi

# Test BIOM
python -c "import biom; print('BIOM:', biom.__version__)" &>> "$LOG_FILE"
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ BIOM-format: INSTALLED${NC}"
else
    echo -e "${YELLOW}! BIOM-format: Check manually${NC}"
fi

# Test Streamlit
python -c "import streamlit; print('Streamlit:', streamlit.__version__)" &>> "$LOG_FILE"
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Streamlit: INSTALLED${NC}"
else
    echo -e "${YELLOW}! Streamlit: Check manually${NC}"
fi

echo -e "\n${BLUE}===============================================${NC}"
echo -e "${GREEN}SETUP COMPLETE!${NC}"
echo -e "\n${YELLOW}Next steps:${NC}"
echo -e "1. Activate: ${GREEN}conda activate $ENV_NAME${NC}"
echo -e "2. Test DADA2: ${GREEN}Rscript -e 'library(dada2)'${NC}"
echo -e "3. Test SourceTracker2: ${GREEN}python -c 'import sourcetracker'${NC}"
echo -e "4. Run CLI: ${GREEN}./run_sourcetracker.sh <fastq_dir> [output_dir] [sources]${NC}"
echo -e "5. Run GUI: ${GREEN}./run_app.sh${NC}"
echo -e "\n${YELLOW}Logs saved to: ${NC}$LOG_FILE"
echo -e "${BLUE}===============================================${NC}"

conda deactivate
