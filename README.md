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

# All source categories (default: V3-V4 amplicons)
./run_sourcetracker.sh <fastq_dir> [output_dir]

# Specific source categories only (available sources: human,cow,pig,duck,chicken,groundwater,river,seawater)
./run_sourcetracker.sh <fastq_dir> [output_dir] "human,cow,pig,duck,chicken,river"

# Specify amplicon type (v4, v34, or v45)
./run_sourcetracker.sh <fastq_dir> [output_dir] "all" v4      # V4 only amplicons
./run_sourcetracker.sh <fastq_dir> [output_dir] "all" v34     # V3-V4 amplicons (default)
./run_sourcetracker.sh <fastq_dir> [output_dir] "all" v45     # V4-V5 amplicons

# Learn error rates from input data (slower but more accurate for different sequencing runs)
./run_sourcetracker.sh <fastq_dir> [output_dir] "all" v34 learn
```

The `fastq_dir` should contain paired FASTQ files named `sample_R1.fastq.gz` / `sample_R2.fastq.gz` (or `sample_1.fastq.gz` / `sample_2.fastq.gz`).

### Input Data Requirements

The pipeline accepts paired-end MiSeq 16S rRNA amplicon data. The reference database uses the V4 region, so V4 is extracted from your input reads.

**Supported amplicon types:**

| Type | Description | R1 Contains | R2 Contains |
|------|-------------|-------------|-------------|
| `v34` | V3-V4 (default) | V3 + 515F + V4 | 806R + V4 |
| `v4` | V4 only | 515F + V4 | 806R + V4 |
| `v45` | V4-V5 | 515F + V4 + V5 | 926R + V5 |

**Primer handling by amplicon type:**

- **V3-V4 (`v34`)**:
  - R1: Searches for 515F anywhere (required to locate V4 start)
  - R2: Trims 806R if found, uses as-is if not found

- **V4 (`v4`)**:
  - R1: Searches for 515F at start, uses as-is if not found (primers may be pre-trimmed)
  - R2: Trims 806R if found, uses as-is if not found

- **V4-V5 (`v45`)**:
  - R1: Searches for 515F at start, uses as-is if not found
  - R2: Contains V5 region (not V4), used as-is with warning
  - ⚠️ Limited support: R2 won't match V4 database well

**Default primers:**
- 515F (V4 forward): `GTGYCAGCMGCCGCGGTAA`
- 806R (V4 reverse): `GACTACNVGGGTWTCTAAT`

This flexibility allows the tool to work with both raw amplicon data (primers intact) and pre-processed data (primers already trimmed).

**DADA2 Error Rates:**

The pipeline offers two options for DADA2 error rate estimation:

| Option | Speed | When to Use |
|--------|-------|-------------|
| Pre-computed (default) | Fast | Data from similar sequencing runs as the reference database |
| Learn from input | Slow | Data from different sequencing platforms or runs |

The pre-computed error rates were generated from V3-V4 amplicon data. For V4-only or V4-V5 data, or data from different sequencing runs, consider using the "learn" option for better accuracy.

### Streamlit GUI

**Windows (double-click):**
- Double-click `Launch_MST_GUI.bat` to start the GUI
- Or run `Create_Desktop_Shortcut.bat` once to add a shortcut to your desktop

**Linux/WSL terminal:**
```bash
./run_app.sh
```

Opens a web interface at http://localhost:8501 where you can:
- Upload FASTQ files
- Select amplicon type (V4, V3-V4, or V4-V5)
- Choose to learn error rates from input data (slower but more accurate)
- Choose source categories
- Visualize results

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
├── setup.sh                    # Environment setup script
├── run_sourcetracker.sh        # CLI pipeline
├── run_app.sh                  # Streamlit GUI launcher (Linux/WSL)
├── Launch_MST_GUI.bat          # Streamlit GUI launcher (Windows double-click)
├── Create_Desktop_Shortcut.bat # Creates desktop shortcut (Windows)
├── app.py                      # Streamlit application
├── requirements.txt            # Python pip dependencies
├── scripts/
│   └── get_v4_from_v34.py      # V4 region extraction
└── db/                         # Reference database
```
