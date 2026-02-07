# SourceTracker2 MST Pipeline

Microbial Source Tracking pipeline combining DADA2 ASV inference with SourceTracker2 Bayesian source estimation. Processes paired-end 16S rRNA amplicon data (V3-V4, V4, or V4-V5 region) and estimates the proportions of microbial communities originating from different source environments.

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

# Basic usage (V3-V4 amplicons, all sources)
./run_sourcetracker.sh <fastq_dir> [output_dir]

# Specify amplicon type
./run_sourcetracker.sh <fastq_dir> <output_dir> all v34    # V3-V4 (default)
./run_sourcetracker.sh <fastq_dir> <output_dir> all v4     # V4 only
./run_sourcetracker.sh <fastq_dir> <output_dir> all v45    # V4-V5

# Specific source categories (available: human,cow,pig,duck,chicken,groundwater,river,seawater)
./run_sourcetracker.sh <fastq_dir> <output_dir> "human,cow,pig"

# Learn error rates from input data (slower but more accurate)
./run_sourcetracker.sh <fastq_dir> <output_dir> all v34 learn
```

The `fastq_dir` should contain paired FASTQ files named `sample_R1.fastq.gz` / `sample_R2.fastq.gz` (or `sample_1.fastq.gz` / `sample_2.fastq.gz`).

### Input Data Requirements

The pipeline accepts paired-end MiSeq data from multiple amplicon types. The reference database uses the V4 region, so V4 is extracted from your input reads.

**Supported amplicon types:**

| Amplicon | Description | R1 Processing | R2 Processing |
|----------|-------------|---------------|---------------|
| **v34** | V3-V4 (~460bp) | Find 515F anywhere, keep after | Find 806R, trim, truncate 150bp |
| **v4** | V4 only (~253bp) | Find 515F at 5' end, trim | Find 806R, trim, truncate 150bp |
| **v45** | V4-V5 (~410bp) | Find 515F, trim; find 806R, truncate before | Find 926R, trim; find 806R_rc, keep after |

**Primers used:**
- 515F (V4 forward): `GTGYCAGCMGCCGCGGTAA` (19bp)
- 806R (V4 reverse): `GACTACNVGGGTWTCTAAT` (19bp)
- 926R (V5 reverse): `CCGYCAATTYMTTTRAGTTT` (20bp)

**Auto-detection:** The pipeline automatically detects whether primers are present or already trimmed:
- If primers found → trims them
- If primers not found → assumes already trimmed, uses reads as-is

This flexibility allows the tool to work with both raw amplicon data (primers intact) and pre-processed data (primers already trimmed).

### Streamlit GUI

```bash
./run_app.sh
```

Opens a web interface at http://localhost:8501 where you can upload FASTQ files, select source categories, and visualize results.

## Pipeline Steps

1. **V4 Extraction** - Extracts V4 region from amplicons using primer search. Supports V3-V4, V4, and V4-V5 input data. Auto-detects primer presence.
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
│   └── get_v4_from_all.py  # V4 region extraction (supports v34, v4, v45)
└── db/                   # Reference database
```
