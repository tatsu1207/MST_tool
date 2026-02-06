# SourceTracker2 MST Pipeline

Microbial Source Tracking pipeline combining DADA2 ASV inference with SourceTracker2 Bayesian source estimation. Processes paired-end 16S rRNA amplicon data (V3-V4 or V4 region) and estimates the proportions of microbial communities originating from different source environments.

## Setup

### Prerequisites
- Linux or WSL2 (Windows Subsystem for Linux)
- [Miniforge](https://github.com/conda-forge/miniforge) (recommended) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

### WSL2 Setup (Windows users)

1. Open PowerShell as Administrator and install WSL2:
   ```powershell
   wsl --install -d Ubuntu
   ```
2. Launch Ubuntu from the Start menu and set up your username/password.
3. Install Miniforge inside WSL:
   ```bash
   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
   bash Miniforge3-Linux-x86_64.sh
   # Restart your terminal after installation
   ```

### Install

```bash
git clone <repo-url>
cd MST_tool_ST
bash setup.sh
```

This creates a conda environment `dada2_mst` with all dependencies (DADA2, SourceTracker2, Streamlit, etc.).

## Usage

### Command Line

```bash
conda activate dada2_mst

# All source categories
./run_sourcetracker.sh <fastq_dir> [output_dir]

# Specific source categories only (available sources:human,cow,pig,duck,chicken,groundwater,river,seawater )
./run_sourcetracker.sh <fastq_dir> [output_dir] "human,cow,pig,duck,chicken,river" 
```

The `fastq_dir` should contain paired FASTQ files named `sample_R1.fastq.gz` / `sample_R2.fastq.gz` (or `sample_1.fastq.gz` / `sample_2.fastq.gz`).

### Input Data Requirements

The pipeline accepts V3-V4 region paired-end MiSeq data. The reference database uses the V4 region, so V4 is extracted from your input reads.

**Primer handling:**
- **R1**: Must contain the V4 forward primer (515F: `GTGYCAGCMGCCGCGGTAA`). The script searches for this primer and keeps everything after it.
- **R2**: Can have primers intact OR already removed:
  - If 806R primer (`GACTACNVGGGTWTCTAAT`) is found → trimmed automatically
  - If 806R primer is not found → assumes primers were already removed and uses R2 as-is
  - R2 is truncated to 150bp in both cases

This flexibility allows the tool to work with both raw amplicon data (primers intact) and pre-processed data (primers already trimmed).

### Streamlit GUI

```bash
./run_app.sh
```

Opens a web interface at http://localhost:8501 where you can upload FASTQ files, select source categories, and visualize results.

## Pipeline Steps

1. **V4 Extraction** - Extracts V4 region from V3-V4 amplicons using primer search (515F/806R). Handles both primer-intact and primer-trimmed input data.
2. **DADA2 Processing** - Filters, denoises, merges reads, and removes chimeras using pre-computed error rates
3. **ASV Mapping** - Maps sink ASVs to the reference database by exact sequence match
4. **Feature Table** - Creates combined source + sink abundance table with rare ASV filtering
5. **SourceTracker2** - Runs Gibbs sampling to estimate source proportions

## Database

The `db/` directory contains the reference source database:

| File | Description |
|------|-------------|
| `db.fasta` | Reference ASV sequences (152K ASVs) |
| `db_table.csv.gz` | ASV abundance matrix (gzipped) |
| `db.design` | Sample-to-environment mapping |
| `err_F.rds` / `err_R.rds` | Pre-computed DADA2 error rates |

## File Structure

```
├── setup.sh              # Environment setup script
├── run_sourcetracker.sh  # CLI pipeline
├── run_app.sh            # Streamlit GUI launcher
├── app.py                # Streamlit application
├── requirements.txt      # Python pip dependencies
├── scripts/
│   └── get_v4_from_v34.py  # V4 region extraction
└── db/                   # Reference database
```
