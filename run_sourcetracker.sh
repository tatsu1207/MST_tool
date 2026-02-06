#!/bin/bash

# SourceTracker2 Pipeline Script
# Usage: ./run_sourcetracker.sh <fastq_dir> [output_dir] [sources]
#
# Arguments:
#   fastq_dir  - Directory containing paired FASTQ files (required)
#                Files should be named: sample_R1.fastq.gz / sample_R2.fastq.gz
#   output_dir - Output directory (optional, default: ./st2_output)
#   sources    - Comma-separated source categories to include (optional, default: all)
#                e.g. "human,cow,pig" or "human,chicken,seawater"

set -e

# Parse arguments
FASTQ_DIR="$1"
OUTPUT_DIR="${2:-./st2_output}"
SOURCES="${3:-all}"

# Validate arguments
if [ -z "$FASTQ_DIR" ]; then
    echo "Usage: $0 <fastq_dir> [output_dir] [sources]"
    echo ""
    echo "Arguments:"
    echo "  fastq_dir  - Directory containing paired FASTQ files (required)"
    echo "               Files should be named: sample_R1.fastq.gz / sample_R2.fastq.gz"
    echo "  output_dir - Output directory (optional, default: ./st2_output)"
    echo "  sources    - Comma-separated source categories (optional, default: all)"
    echo "               e.g. \"human,cow,pig\""
    exit 1
fi

if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

# Find and pair R1/R2 files using Python for robust pattern matching
eval "$(python3 - "$FASTQ_DIR" << 'PAIR_SCRIPT'
import sys, re, os
fastq_dir = sys.argv[1]
files = sorted(f for f in os.listdir(fastq_dir) if f.endswith(('.fastq.gz', '.fq.gz')))
pattern = re.compile(r'^(.+?)[._]([Rr]?[12])(?:[._](?:001|002))?\.(?:fastq|fq)\.gz$')
pairs = {}
for f in files:
    m = pattern.match(f)
    if not m:
        continue
    base, tag = m.group(1), m.group(2)
    read = 'R1' if tag in ('R1','r1','1') else 'R2' if tag in ('R2','r2','2') else None
    if read:
        pairs.setdefault(base, {})
        pairs[base][read] = f
# Output bash arrays for complete pairs only
r1_list, r2_list, name_list = [], [], []
for base, p in sorted(pairs.items()):
    if 'R1' in p and 'R2' in p:
        r1_list.append(os.path.join(fastq_dir, p['R1']))
        r2_list.append(os.path.join(fastq_dir, p['R2']))
        name_list.append(base)
print(f'R1_FILES=({" ".join(repr(f) for f in r1_list)})')
print(f'R2_FILES=({" ".join(repr(f) for f in r2_list)})')
print(f'SAMPLE_NAMES=({" ".join(repr(n) for n in name_list)})')
PAIR_SCRIPT
)"

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "Error: No paired FASTQ files found in $FASTQ_DIR"
    echo "Expected naming: sample_R1.fastq.gz / sample_R2.fastq.gz"
    echo "            or: sample_1.fastq.gz / sample_2.fastq.gz"
    exit 1
fi

# Get script directory (where db files are located)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_FASTA="$SCRIPT_DIR/db/db.fasta"
DB_TABLE="$SCRIPT_DIR/db/db_table.csv.gz"
DB_DESIGN="$SCRIPT_DIR/db/db.design"

# Pre-computed error rate files
ERR_F_FILE="$SCRIPT_DIR/db/err_F.rds"
ERR_R_FILE="$SCRIPT_DIR/db/err_R.rds"

# Validate database files
for f in "$DB_FASTA" "$DB_TABLE" "$DB_DESIGN"; do
    if [ ! -f "$f" ]; then
        echo "Error: Database file not found: $f"
        exit 1
    fi
done

if [ ! -f "$ERR_F_FILE" ] || [ ! -f "$ERR_R_FILE" ]; then
    echo "Error: Pre-computed error rate files not found"
    echo "Expected: $ERR_F_FILE and $ERR_R_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "SourceTracker2 Pipeline"
echo "========================================"
echo "FASTQ directory: $FASTQ_DIR"
echo "Samples found: ${#SAMPLE_NAMES[@]}"
for S in "${SAMPLE_NAMES[@]}"; do
    echo "  - $S"
done
echo "Output directory: $OUTPUT_DIR"
echo "Sources: $SOURCES"
echo "Error rates: pre-computed"
echo "========================================"

# Activate conda environment
echo ""
echo "[Step 1/6] Activating conda environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dada2_mst

# V4 extraction script
V4_EXTRACT_SCRIPT="$SCRIPT_DIR/scripts/get_v4_from_v34.py"

TRIMMED_DIR="$(realpath "$OUTPUT_DIR")/trimmed"
mkdir -p "$TRIMMED_DIR"
WORK_DIR="$(realpath "$OUTPUT_DIR")/work"
mkdir -p "$WORK_DIR"

# Process each sample through V4 extraction and DADA2
SINK_ABUNDANCE_FILES=()
PROCESSED_NAMES=()
TOTAL=${#R1_FILES[@]}

for IDX in "${!R1_FILES[@]}"; do
    R1="${R1_FILES[$IDX]}"
    R2="${R2_FILES[$IDX]}"
    SAMPLE_NAME="${SAMPLE_NAMES[$IDX]}"
    NUM=$((IDX + 1))

    echo ""
    echo "========================================"
    echo "Processing sample $NUM/$TOTAL: $SAMPLE_NAME"
    echo "========================================"

    # Step 2: Extract V4 region
    echo ""
    echo "[Step 2/6] Extracting V4 region from $SAMPLE_NAME..."

    R1_ABS=$(realpath "$R1")
    R2_ABS=$(realpath "$R2")
    # Use symlinks if possible, fall back to copy (WSL may not support symlinks)
    ln -sf "$R1_ABS" "$WORK_DIR/${SAMPLE_NAME}_R1.fastq.gz" 2>/dev/null || \
        cp -f "$R1_ABS" "$WORK_DIR/${SAMPLE_NAME}_R1.fastq.gz"
    ln -sf "$R2_ABS" "$WORK_DIR/${SAMPLE_NAME}_R2.fastq.gz" 2>/dev/null || \
        cp -f "$R2_ABS" "$WORK_DIR/${SAMPLE_NAME}_R2.fastq.gz"

    cd "$WORK_DIR"
    python3 "$V4_EXTRACT_SCRIPT" "$SAMPLE_NAME" "$TRIMMED_DIR" --mismatches 2
    cd - > /dev/null

    TRIMMED_R1="$TRIMMED_DIR/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
    TRIMMED_R2="$TRIMMED_DIR/${SAMPLE_NAME}_trimmed_R2.fastq.gz"

    if [ ! -f "$TRIMMED_R1" ] || [ ! -f "$TRIMMED_R2" ]; then
        echo "Warning: V4 extraction failed for $SAMPLE_NAME, skipping..."
        continue
    fi

    echo "V4 extraction complete."

    # Step 3: Process with DADA2
    echo ""
    echo "[Step 3/6] Processing $SAMPLE_NAME with DADA2..."

    Rscript --vanilla - "$TRIMMED_R1" "$TRIMMED_R2" "$OUTPUT_DIR" "$SAMPLE_NAME" "$DB_FASTA" "$ERR_F_FILE" "$ERR_R_FILE" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
r1_file <- args[1]
r2_file <- args[2]
output_dir <- args[3]
sample_name <- args[4]
db_fasta <- args[5]
err_F_file <- args[6]
err_R_file <- args[7]

library(dada2)
library(Biostrings)

cat("Processing paired-end reads (V4 region extracted)...\n")

# Create temp directory for processing
temp_dir <- file.path(output_dir, paste0("temp_dada2_", sample_name))
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

# Filter and trim sink sample
filt_r1 <- file.path(temp_dir, "filtered_R1.fastq.gz")
filt_r2 <- file.path(temp_dir, "filtered_R2.fastq.gz")

cat("Filtering reads (V4 already extracted by primer search)...\n")

filt_out <- filterAndTrim(
    r1_file, filt_r1,
    r2_file, filt_r2,
    truncLen = c(0, 0),
    maxN = 0,
    maxEE = c(2, 2),
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
)

cat(sprintf("Input reads: %d, Filtered reads: %d\n", filt_out[1,1], filt_out[1,2]))

if (filt_out[1,2] == 0) {
    stop("No reads passed the filter!")
}

# Load pre-computed error rates
cat("Loading pre-computed error rates...\n")
err_F <- readRDS(err_F_file)
err_R <- readRDS(err_R_file)

# Dereplicate sink sample
cat("Dereplicating sink sample...\n")
derep_r1 <- derepFastq(filt_r1, verbose = FALSE)
derep_r2 <- derepFastq(filt_r2, verbose = FALSE)

# Denoise using database error rates
cat("Denoising with database error rates...\n")
dada_r1 <- dada(derep_r1, err = err_F, multithread = TRUE, verbose = FALSE)
dada_r2 <- dada(derep_r2, err = err_R, multithread = TRUE, verbose = FALSE)

# Merge paired reads
cat("Merging paired reads...\n")
merged <- mergePairs(dada_r1, derep_r1, dada_r2, derep_r2, verbose = TRUE)

# Check merge results
if (is.data.frame(merged)) {
    n_merged <- nrow(merged[merged$accept, ])
} else {
    n_merged <- sum(sapply(merged, function(x) nrow(x[x$accept, ])))
}
cat(sprintf("Successfully merged reads: %d\n", n_merged))

if (n_merged == 0) {
    stop("No reads were successfully merged! Check if reads overlap properly.")
}

# Make sequence table
seqtab <- makeSequenceTable(merged)
cat(sprintf("Sequence table dimensions: %d samples x %d ASVs\n", nrow(seqtab), ncol(seqtab)))

if (ncol(seqtab) == 0) {
    stop("No ASVs in sequence table after merging!")
}

# Remove chimeras
cat("Removing chimeras...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = FALSE)
cat(sprintf("ASVs after chimera removal: %d\n", ncol(seqtab_nochim)))

if (ncol(seqtab_nochim) == 0) {
    stop("No ASVs remained after chimera removal!")
}

# Get sink sequences
sink_seqs <- colnames(seqtab_nochim)
sink_counts <- as.vector(seqtab_nochim[1,])

# Read database sequences
cat("Loading database sequences...\n")
db_seqs <- readDNAStringSet(db_fasta)
db_names <- names(db_seqs)
db_sequences <- as.character(db_seqs)

# Map sink ASVs to database ASVs by exact sequence match
cat("Mapping sink ASVs to database...\n")
sink_to_db <- data.frame(
    sink_seq = sink_seqs,
    sink_count = sink_counts,
    db_asv = NA_character_,
    stringsAsFactors = FALSE
)

for (i in seq_along(sink_seqs)) {
    match_idx <- which(db_sequences == sink_seqs[i])
    if (length(match_idx) > 0) {
        sink_to_db$db_asv[i] <- db_names[match_idx[1]]
    }
}

# Summary of mapping
mapped_count <- sum(!is.na(sink_to_db$db_asv))
mapped_reads <- sum(sink_to_db$sink_count[!is.na(sink_to_db$db_asv)])
total_reads <- sum(sink_to_db$sink_count)
cat(sprintf("Mapped ASVs: %d/%d (%.1f%%)\n", mapped_count, nrow(sink_to_db), 100*mapped_count/nrow(sink_to_db)))
cat(sprintf("Mapped reads: %d/%d (%.1f%%)\n", mapped_reads, total_reads, 100*mapped_reads/total_reads))

# Create sink abundance vector for database ASVs
sink_abundance <- setNames(rep(0, length(db_names)), db_names)
for (i in seq_len(nrow(sink_to_db))) {
    if (!is.na(sink_to_db$db_asv[i])) {
        asv_name <- sink_to_db$db_asv[i]
        sink_abundance[asv_name] <- sink_abundance[asv_name] + sink_to_db$sink_count[i]
    }
}

# Save sink abundance
sink_file <- file.path(output_dir, paste0(sample_name, "_abundance.tsv"))
write.table(
    data.frame(ASV = names(sink_abundance), abundance = as.numeric(sink_abundance)),
    file = sink_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

# Save mapping summary
mapping_file <- file.path(output_dir, paste0(sample_name, "_mapping_summary.tsv"))
write.table(sink_to_db, file = mapping_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("Sink abundance saved to: %s\n", sink_file))
cat(sprintf("Mapping summary saved to: %s\n", mapping_file))

# Clean up temp directory
unlink(temp_dir, recursive = TRUE)

cat("DADA2 processing complete.\n")
RSCRIPT

    SINK_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_abundance.tsv"
    if [ -f "$SINK_FILE" ]; then
        SINK_ABUNDANCE_FILES+=("$SINK_FILE")
        PROCESSED_NAMES+=("$SAMPLE_NAME")
    else
        echo "Warning: DADA2 processing failed for $SAMPLE_NAME, skipping..."
    fi
done

# Check that at least one sample succeeded
if [ ${#PROCESSED_NAMES[@]} -eq 0 ]; then
    echo "Error: No samples were successfully processed."
    exit 1
fi

echo ""
echo "========================================"
echo "Successfully processed ${#PROCESSED_NAMES[@]}/${TOTAL} samples"
echo "========================================"

# Build comma-separated lists for python scripts
SINK_FILES_CSV=$(IFS=,; echo "${SINK_ABUNDANCE_FILES[*]}")
SINK_NAMES_CSV=$(IFS=,; echo "${PROCESSED_NAMES[*]}")

# Step 4: Create combined feature table for SourceTracker2
echo ""
echo "[Step 4/6] Creating combined feature table..."

python3 - "$DB_TABLE" "$OUTPUT_DIR" "$DB_DESIGN" "$SOURCES" "$SINK_FILES_CSV" "$SINK_NAMES_CSV" << 'PYTHON_SCRIPT'
import sys
import pandas as pd

db_table_file = sys.argv[1]
output_dir = sys.argv[2]
db_design_file = sys.argv[3]
sources_arg = sys.argv[4]
sink_files = sys.argv[5].split(',')
sink_names = sys.argv[6].split(',')

print("Loading database table...")
# Read database table (ASVs as rows, samples as columns)
db_table = pd.read_csv(db_table_file, index_col=0)
print(f"Database table: {db_table.shape[0]} ASVs x {db_table.shape[1]} samples")

# Filter by selected source categories
design = pd.read_csv(db_design_file, sep='\t')
design.columns = ['SampleID', 'Env']
if sources_arg != "all":
    selected_sources = [s.strip() for s in sources_arg.split(',')]
    selected_samples = design[design['Env'].isin(selected_sources)]['SampleID'].tolist()
else:
    selected_samples = design['SampleID'].tolist()
available_cols = [col for col in selected_samples if col in db_table.columns]
db_table = db_table[available_cols]
print(f"After source filtering: {db_table.shape[1]} samples")

# Add all sink samples
combined = db_table.copy()
for sink_file, sample_name in zip(sink_files, sink_names):
    print(f"Loading sink abundance for {sample_name}...")
    sink_df = pd.read_csv(sink_file, sep='\t')
    sink_abundance = dict(zip(sink_df['ASV'], sink_df['abundance']))
    sink_col = pd.Series([sink_abundance.get(asv, 0) for asv in combined.index],
                          index=combined.index, name=sample_name)
    combined[sample_name] = sink_col

# Filter rare ASVs to reduce matrix size and speed up SourceTracker2
n_before = len(combined)
min_samples = 2
min_total = 10
mask = (combined > 0).sum(axis=1) >= min_samples
mask &= combined.sum(axis=1) >= min_total
combined = combined[mask]
print(f"ASV filtering: {n_before} -> {len(combined)} (removed {n_before - len(combined)} rare ASVs)")

# Save as BIOM-compatible TSV (ASVs as rows, samples as columns, #OTU ID index)
combined.index.name = '#OTU ID'

combined_file = f"{output_dir}/combined_table.tsv"
combined.to_csv(combined_file, sep='\t')
print(f"Combined table saved: {combined.shape[0]} features x {combined.shape[1]} samples")
print(f"Saved to: {combined_file}")
PYTHON_SCRIPT

# Step 5: Create mapping file for SourceTracker2
echo ""
echo "[Step 5/6] Creating mapping file..."

python3 - "$DB_DESIGN" "$OUTPUT_DIR" "$SOURCES" "$SINK_NAMES_CSV" << 'PYTHON_SCRIPT2'
import sys
import pandas as pd

db_design_file = sys.argv[1]
output_dir = sys.argv[2]
sources_arg = sys.argv[3]
sink_names = sys.argv[4].split(',')

print("Creating SourceTracker2 mapping file...")

# Read design file
design = pd.read_csv(db_design_file, sep='\t')
design.columns = ['SampleID', 'Env']

# Filter by selected source categories
if sources_arg != "all":
    selected_sources = [s.strip() for s in sources_arg.split(',')]
    design = design[design['Env'].isin(selected_sources)]
    print(f"Filtered to sources: {selected_sources}")

# Add SourceSink column - all existing samples are sources
design['SourceSink'] = 'source'

# Add all sink samples
sink_rows = pd.DataFrame({
    'SampleID': sink_names,
    'Env': ['sink'] * len(sink_names),
    'SourceSink': ['sink'] * len(sink_names)
})

mapping = pd.concat([design, sink_rows], ignore_index=True)

# Set index
mapping = mapping.set_index('SampleID')

# Save mapping file
mapping_file = f"{output_dir}/mapping.tsv"
mapping.to_csv(mapping_file, sep='\t')
print(f"Mapping file saved with {len(mapping)} samples ({len(sink_names)} sinks)")
print(f"Source categories: {design['Env'].nunique()}")
print(f"Saved to: {mapping_file}")
PYTHON_SCRIPT2

# Step 6: Run SourceTracker2
echo ""
echo "[Step 6/6] Running SourceTracker2..."

ST2_OUTPUT="$OUTPUT_DIR/sourcetracker2_results"
rm -rf "$ST2_OUTPUT"

sourcetracker2 \
    -i "$OUTPUT_DIR/combined_table.tsv" \
    -m "$OUTPUT_DIR/mapping.tsv" \
    -o "$ST2_OUTPUT" \
    --source_sink_column SourceSink \
    --source_column_value source \
    --sink_column_value sink \
    --source_category_column Env \
    --jobs 4

echo ""
echo "========================================"
echo "Pipeline Complete!"
echo "========================================"
echo ""
echo "Samples processed: ${#PROCESSED_NAMES[@]}"
for S in "${PROCESSED_NAMES[@]}"; do
    echo "  - $S"
done
echo ""
echo "Results saved to: $ST2_OUTPUT"
echo ""
echo "Key output files:"
echo "  - mixing_proportions.txt: Source contributions for sink samples"
echo "  - mixing_proportions_stds.txt: Standard deviations"
echo ""

# Display results
if [ -f "$ST2_OUTPUT/mixing_proportions.txt" ]; then
    echo "Source Contributions:"
    echo "----------------------------------------"
    cat "$ST2_OUTPUT/mixing_proportions.txt"
fi
