#!/usr/bin/env python3
"""
SourceTracker2 Streamlit GUI
Allows users to upload fastq files and select source categories for microbial source tracking.
"""

import streamlit as st
import pandas as pd
import subprocess
import tempfile
import os
import shutil
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

# Configuration
SCRIPT_DIR = Path(__file__).parent.resolve()
DB_FASTA = SCRIPT_DIR / "db" / "db.fasta"
DB_TABLE = SCRIPT_DIR / "db" / "db_table.csv.gz"
DB_DESIGN = SCRIPT_DIR / "db" / "db.design"

# Pre-computed error rate files
ERR_F_FILE = SCRIPT_DIR / "db" / "err_F.rds"
ERR_R_FILE = SCRIPT_DIR / "db" / "err_R.rds"

# Conda environment name
CONDA_ENV = "dada2_mst"

def get_source_categories():
    """Auto-detect source categories from db.design."""
    design = pd.read_csv(DB_DESIGN, sep='\t', header=0, names=['SampleID', 'Group'])
    return sorted(design['Group'].unique().tolist())

SOURCE_CATEGORIES = get_source_categories()

st.set_page_config(
    page_title="SourceTracker2 GUI",
    page_icon="🧬",
    layout="wide"
)

def check_error_rates_exist():
    """Check if pre-computed error rates exist."""
    return ERR_F_FILE.exists() and ERR_R_FILE.exists()

def load_design_file():
    """Load and parse the design file."""
    design = pd.read_csv(DB_DESIGN, sep='\t', header=0, names=['SampleID', 'Group'])
    return design

def get_source_counts(design):
    """Get count of samples per source category."""
    return design['Group'].value_counts().to_dict()

def save_uploaded_files(uploaded_files, temp_dir):
    """Save uploaded files to temporary directory."""
    saved_paths = []
    for uploaded_file in uploaded_files:
        file_path = Path(temp_dir) / uploaded_file.name
        with open(file_path, 'wb') as f:
            f.write(uploaded_file.getbuffer())
        saved_paths.append(file_path)
    return saved_paths

def pair_fastq_files(file_paths):
    """Pair R1 and R2 fastq files based on naming convention."""
    import re
    pairs = {}
    # Patterns anchored before the file extension to avoid matching mid-name
    # Each tuple: (regex, read_number_group_index)
    # Matches _R1/_R2 or _1/_2 only at the end before .fastq.gz/.fq.gz
    pattern = re.compile(r'^(.+?)[._]([Rr]?[12])(?:[._](?:001|002))?\.(?:fastq|fq)\.gz$')

    for fp in file_paths:
        name = fp.name
        m = pattern.match(name)
        if not m:
            continue

        base_name = m.group(1)
        read_tag = m.group(2)

        if read_tag in ('R1', 'r1', '1'):
            read_type = 'R1'
        elif read_tag in ('R2', 'r2', '2'):
            read_type = 'R2'
        else:
            continue

        if base_name not in pairs:
            pairs[base_name] = {'R1': None, 'R2': None, 'sample_name': base_name}
        pairs[base_name][read_type] = fp

    # Filter complete pairs
    complete_pairs = {k: v for k, v in pairs.items() if v['R1'] and v['R2']}
    return complete_pairs

def create_filtered_design(design, selected_sources):
    """Create a filtered design file with only selected sources."""
    filtered = design[design['Group'].isin(selected_sources)]
    return filtered

# V4 extraction script path
V4_EXTRACT_SCRIPT = SCRIPT_DIR / "scripts" / "get_v4_from_all.py"

def run_v4_extraction(r1_path, r2_path, output_dir, sample_name, amplicon="v34"):
    """Extract V4 region from amplicons using primer search."""
    trimmed_dir = Path(output_dir) / "trimmed"
    trimmed_dir.mkdir(exist_ok=True)

    work_dir = Path(output_dir) / "work"
    work_dir.mkdir(exist_ok=True)

    # Create links with expected names for the script
    # Use symlinks if possible, fall back to copy (WSL may not support symlinks)
    r1_link = work_dir / f"{sample_name}_R1.fastq.gz"
    r2_link = work_dir / f"{sample_name}_R2.fastq.gz"

    r1_link.unlink(missing_ok=True)
    r2_link.unlink(missing_ok=True)

    try:
        r1_link.symlink_to(Path(r1_path).resolve())
        r2_link.symlink_to(Path(r2_path).resolve())
    except OSError:
        import shutil
        shutil.copy2(Path(r1_path).resolve(), r1_link)
        shutil.copy2(Path(r2_path).resolve(), r2_link)

    # Run the V4 extraction script
    cmd = f"cd {work_dir} && python3 {V4_EXTRACT_SCRIPT} {sample_name} {trimmed_dir} --mismatches 2 --amplicon {amplicon}"

    result = subprocess.run(
        f"source $(conda info --base)/etc/profile.d/conda.sh && conda activate dada2_mst && {cmd}",
        shell=True, capture_output=True, text=True, executable='/bin/bash'
    )

    if result.returncode != 0:
        return None, None, result

    # Output files from the script
    trimmed_r1 = trimmed_dir / f"{sample_name}_trimmed_R1.fastq.gz"
    trimmed_r2 = trimmed_dir / f"{sample_name}_trimmed_R2.fastq.gz"

    if not trimmed_r1.exists() or not trimmed_r2.exists():
        return None, None, result

    return trimmed_r1, trimmed_r2, result

def run_dada2_processing(r1_path, r2_path, output_dir, sample_name, learn_error_rates=False):
    """Run DADA2 processing on sink sample."""
    # Generate error rate loading/learning code based on option
    if learn_error_rates:
        error_rate_code = '''
# Learn error rates from input data
cat("Learning error rates from input data (this may take a while)...\\n")
err_F <- learnErrors(filt_r1, multithread = TRUE, verbose = TRUE)
err_R <- learnErrors(filt_r2, multithread = TRUE, verbose = TRUE)
cat("Error rate learning complete.\\n")
'''
    else:
        error_rate_code = f'''
# Load pre-computed error rates
cat("Loading pre-computed error rates...\\n")
err_F <- readRDS("{ERR_F_FILE}")
err_R <- readRDS("{ERR_R_FILE}")
'''

    r_script = f'''
library(dada2)
library(Biostrings)

r1_file <- "{r1_path}"
r2_file <- "{r2_path}"
output_dir <- "{output_dir}"
sample_name <- "{sample_name}"
db_fasta <- "{DB_FASTA}"

cat("Processing paired-end reads...\\n")

temp_dir <- file.path(output_dir, "temp_dada2")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

filt_r1 <- file.path(temp_dir, "filtered_R1.fastq.gz")
filt_r2 <- file.path(temp_dir, "filtered_R2.fastq.gz")

cat("Filtering reads (V4 already extracted by primer search)...\\n")

# Use truncLen=0 to keep full length (already trimmed by cutadapt)
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

cat(sprintf("Input reads: %d, Filtered reads: %d\\n", filt_out[1,1], filt_out[1,2]))

if (filt_out[1,2] == 0) {{
    stop("No reads passed the filter!")
}}

{error_rate_code}

cat("Dereplicating...\\n")
derep_r1 <- derepFastq(filt_r1, verbose = FALSE)
derep_r2 <- derepFastq(filt_r2, verbose = FALSE)

cat("Denoising...\\n")
dada_r1 <- dada(derep_r1, err = err_F, multithread = TRUE, verbose = FALSE)
dada_r2 <- dada(derep_r2, err = err_R, multithread = TRUE, verbose = FALSE)

cat("Merging paired reads...\\n")
merged <- mergePairs(dada_r1, derep_r1, dada_r2, derep_r2, verbose = TRUE)

# Check merge results
if (is.data.frame(merged)) {{
    n_merged <- nrow(merged[merged$accept, ])
}} else {{
    n_merged <- sum(sapply(merged, function(x) nrow(x[x$accept, ])))
}}
cat(sprintf("Successfully merged reads: %d\\n", n_merged))

if (n_merged == 0) {{
    stop("No reads were successfully merged! Check if reads overlap properly.")
}}

seqtab <- makeSequenceTable(merged)
cat(sprintf("Sequence table: %d samples x %d ASVs\\n", nrow(seqtab), ncol(seqtab)))

if (ncol(seqtab) == 0) {{
    stop("No ASVs in sequence table after merging!")
}}

cat("Removing chimeras...\\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = FALSE)
cat(sprintf("ASVs after chimera removal: %d\\n", ncol(seqtab_nochim)))

if (ncol(seqtab_nochim) == 0) {{
    stop("No ASVs remained after chimera removal!")
}}

sink_seqs <- colnames(seqtab_nochim)
sink_counts <- as.vector(seqtab_nochim[1,])

cat("Loading database sequences...\\n")
db_seqs <- readDNAStringSet(db_fasta)
db_names <- names(db_seqs)
db_sequences <- as.character(db_seqs)

cat("Mapping sink ASVs to database...\\n")
sink_to_db <- data.frame(
    sink_seq = sink_seqs,
    sink_count = sink_counts,
    db_asv = NA_character_,
    stringsAsFactors = FALSE
)

for (i in seq_along(sink_seqs)) {{
    match_idx <- which(db_sequences == sink_seqs[i])
    if (length(match_idx) > 0) {{
        sink_to_db$db_asv[i] <- db_names[match_idx[1]]
    }}
}}

mapped_count <- sum(!is.na(sink_to_db$db_asv))
mapped_reads <- sum(sink_to_db$sink_count[!is.na(sink_to_db$db_asv)])
total_reads <- sum(sink_to_db$sink_count)
cat(sprintf("Mapped ASVs: %d/%d (%.1f%%)\\n", mapped_count, nrow(sink_to_db), 100*mapped_count/nrow(sink_to_db)))
cat(sprintf("Mapped reads: %d/%d (%.1f%%)\\n", mapped_reads, total_reads, 100*mapped_reads/total_reads))

sink_abundance <- setNames(rep(0, length(db_names)), db_names)
for (i in seq_len(nrow(sink_to_db))) {{
    if (!is.na(sink_to_db$db_asv[i])) {{
        asv_name <- sink_to_db$db_asv[i]
        sink_abundance[asv_name] <- sink_abundance[asv_name] + sink_to_db$sink_count[i]
    }}
}}

sink_file <- file.path(output_dir, paste0(sample_name, "_abundance.tsv"))
write.table(
    data.frame(ASV = names(sink_abundance), abundance = as.numeric(sink_abundance)),
    file = sink_file,
    sep = "\\t",
    row.names = FALSE,
    quote = FALSE
)

mapping_file <- file.path(output_dir, paste0(sample_name, "_mapping_summary.tsv"))
write.table(sink_to_db, file = mapping_file, sep = "\\t", row.names = FALSE, quote = FALSE)

unlink(temp_dir, recursive = TRUE)

cat("DADA2 processing complete.\\n")

# Return mapping stats
cat(sprintf("MAPPING_STATS:%d,%d,%d,%d\\n", mapped_count, nrow(sink_to_db), mapped_reads, total_reads))
'''

    # Write R script to temp file
    r_script_path = Path(output_dir) / "process_sink.R"
    with open(r_script_path, 'w') as f:
        f.write(r_script)

    # Run R script
    cmd = f"source $(conda info --base)/etc/profile.d/conda.sh && conda activate dada2_mst && Rscript {r_script_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable='/bin/bash')

    return result

def create_combined_table(sink_files, output_dir, selected_sources):
    """Create combined feature table for SourceTracker2."""
    # Load database table
    db_table = pd.read_csv(DB_TABLE, index_col=0)

    # Load design and filter by selected sources
    design = load_design_file()
    selected_samples = design[design['Group'].isin(selected_sources)]['SampleID'].tolist()

    # Filter database table to selected source samples
    available_cols = [col for col in selected_samples if col in db_table.columns]
    db_filtered = db_table[available_cols]

    # Add sink samples
    combined = db_filtered.copy()
    for sink_file, sample_name in sink_files:
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

    combined_file = Path(output_dir) / "combined_table.tsv"
    combined.to_csv(combined_file, sep='\t')

    return combined_file, len(available_cols)

def create_mapping_file(sink_sample_names, output_dir, selected_sources):
    """Create mapping file for SourceTracker2."""
    design = load_design_file()

    # Filter to selected sources
    filtered_design = design[design['Group'].isin(selected_sources)].copy()
    filtered_design.columns = ['SampleID', 'Env']
    filtered_design['SourceSink'] = 'source'

    # Add sink samples
    sink_rows = pd.DataFrame({
        'SampleID': sink_sample_names,
        'Env': ['sink'] * len(sink_sample_names),
        'SourceSink': ['sink'] * len(sink_sample_names)
    })

    mapping = pd.concat([filtered_design, sink_rows], ignore_index=True)
    mapping = mapping.set_index('SampleID')

    mapping_file = Path(output_dir) / "mapping.tsv"
    mapping.to_csv(mapping_file, sep='\t')

    return mapping_file

def run_sourcetracker2(combined_table, mapping_file, output_dir):
    """Run SourceTracker2 gibbs sampling."""
    st2_output = Path(output_dir) / "sourcetracker2_results"
    # Don't pre-create: sourcetracker2 CLI expects to create this dir itself
    if st2_output.exists():
        import shutil
        shutil.rmtree(st2_output)

    cmd = f"""source $(conda info --base)/etc/profile.d/conda.sh && \
conda activate dada2_mst && \
sourcetracker2 \
    -i {combined_table} \
    -m {mapping_file} \
    -o {st2_output} \
    --source_sink_column SourceSink \
    --source_column_value source \
    --sink_column_value sink \
    --source_category_column Env \
    --jobs 4"""

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable='/bin/bash')

    return st2_output, result

def load_results(results_dir):
    """Load SourceTracker2 results."""
    proportions_file = results_dir / "mixing_proportions.txt"
    stds_file = results_dir / "mixing_proportions_stds.txt"

    if proportions_file.exists():
        proportions = pd.read_csv(proportions_file, sep='\t', index_col=0)
        stds = pd.read_csv(stds_file, sep='\t', index_col=0) if stds_file.exists() else None
        # Ensure samples are rows and sources are columns.
        # Source categories are known, so if they appear as index, transpose.
        source_cats = set(SOURCE_CATEGORIES) | {'Unknown'}
        if source_cats & set(proportions.index):
            proportions = proportions.T
            if stds is not None:
                stds = stds.T
        return proportions, stds
    return None, None

def plot_source_contributions(proportions, sample_name=None):
    """Create a bar plot of source contributions."""
    if sample_name:
        data = proportions.loc[sample_name]
    else:
        data = proportions.iloc[0]

    # Remove 'Unknown' for cleaner visualization if it's very small
    df = pd.DataFrame({
        'Source': data.index,
        'Proportion': data.values
    })
    df = df.sort_values('Proportion', ascending=True)

    fig = px.bar(
        df,
        x='Proportion',
        y='Source',
        orientation='h',
        title=f'Source Contributions{" - " + sample_name if sample_name else ""}',
        color='Proportion',
        color_continuous_scale='Viridis'
    )
    fig.update_layout(
        xaxis_title='Proportion',
        yaxis_title='Source',
        showlegend=False,
        height=400
    )
    return fig

def plot_all_samples_heatmap(proportions):
    """Create a heatmap of all samples."""
    fig = px.imshow(
        proportions,
        labels=dict(x="Source", y="Sample", color="Proportion"),
        title="Source Contributions Heatmap",
        color_continuous_scale='Viridis',
        aspect='auto'
    )
    fig.update_layout(height=max(400, len(proportions) * 30))
    return fig

def plot_pie_chart(proportions, sample_name=None):
    """Create a pie chart of source contributions."""
    if sample_name:
        data = proportions.loc[sample_name]
    else:
        data = proportions.iloc[0]

    # Filter out very small contributions for cleaner pie
    df = pd.DataFrame({
        'Source': data.index,
        'Proportion': data.values
    })
    df = df[df['Proportion'] > 0.001]  # Filter out < 0.1%

    fig = px.pie(
        df,
        values='Proportion',
        names='Source',
        title=f'Source Contributions{" - " + sample_name if sample_name else ""}'
    )
    return fig

# Main App
def main():
    st.title("🧬 SourceTracker2 GUI")
    st.markdown("**Microbial Source Tracking Analysis**")

    # Check error rates exist
    if not check_error_rates_exist():
        st.error("⚠️ Error rate files not found. Please ensure err_F.rds and err_R.rds exist in the script directory.")
        st.stop()

    # Sidebar for configuration
    with st.sidebar:
        st.header("⚙️ Configuration")

        # Amplicon type selection
        st.subheader("Amplicon Region")
        amplicon_options = {
            "v34": "V3-V4 (default)",
            "v4": "V4 only",
            "v45": "V4-V5"
        }
        selected_amplicon = st.selectbox(
            "Select amplicon type",
            options=list(amplicon_options.keys()),
            format_func=lambda x: amplicon_options[x],
            index=0,
            help="V3-V4: primers in middle of R1. V4: primers at start. V4-V5: merged reads are ~410bp (limited support)."
        )

        if selected_amplicon == "v45":
            st.info("ℹ️ V4-V5: V4 region will be extracted by finding 806R boundary (removes V5).")

        st.divider()

        # Error rate option
        st.subheader("DADA2 Error Rates")
        learn_error_rates = st.checkbox(
            "Learn error rates from input data",
            value=False,
            help="If checked, DADA2 will learn error rates from your input files instead of using pre-computed rates. More accurate but significantly slower."
        )
        if learn_error_rates:
            st.info("ℹ️ Learning error rates adds significant processing time but may improve accuracy, especially for data from different sequencing runs.")
        else:
            st.caption("Using pre-computed error rates (fast)")

        st.divider()

        # Source selection
        st.subheader("Select Source Categories")
        design = load_design_file()
        source_counts = get_source_counts(design)

        selected_sources = []
        for source in SOURCE_CATEGORIES:
            count = source_counts.get(source, 0)
            if st.checkbox(f"{source.capitalize()} ({count} samples)", value=True, key=f"src_{source}"):
                selected_sources.append(source)

        if not selected_sources:
            st.warning("Please select at least one source category.")

        st.divider()

        # Database info
        st.subheader("📊 Database Info")
        st.write(f"Total samples: {len(design)}")
        st.write(f"Source categories: {len(SOURCE_CATEGORIES)}")

        # Check database files
        db_status = all([DB_FASTA.exists(), DB_TABLE.exists(), DB_DESIGN.exists()])
        if db_status:
            st.success("✅ Database files found")
        else:
            st.error("❌ Database files missing")


    # Main content
    tab1, tab2 = st.tabs(["📤 Upload & Run", "ℹ️ Help"])

    with tab1:
        st.header("Upload FASTQ Files")
        st.markdown("""
        Upload paired-end FASTQ files (R1 and R2) for your sink samples.
        Files should follow naming convention: `samplename_R1.fastq.gz` and `samplename_R2.fastq.gz`
        """)

        uploaded_files = st.file_uploader(
            "Choose FASTQ files",
            type=['gz', 'fastq'],
            accept_multiple_files=True,
            help="Upload paired-end FASTQ files (.fastq.gz)"
        )

        if uploaded_files:
            # Create temporary directory for processing
            if 'temp_dir' not in st.session_state:
                st.session_state.temp_dir = tempfile.mkdtemp()

            # Save uploaded files
            saved_paths = save_uploaded_files(uploaded_files, st.session_state.temp_dir)

            # Pair files
            pairs = pair_fastq_files(saved_paths)

            if pairs:
                st.success(f"Found {len(pairs)} paired sample(s)")

                # Display paired files
                pair_data = []
                for base, info in pairs.items():
                    pair_data.append({
                        'Sample': info['sample_name'],
                        'R1': info['R1'].name if info['R1'] else 'Missing',
                        'R2': info['R2'].name if info['R2'] else 'Missing'
                    })
                st.dataframe(pd.DataFrame(pair_data), width='stretch')

                # Run analysis button
                if st.button("🚀 Run SourceTracker2", type="primary", disabled=not selected_sources):
                    if not selected_sources:
                        st.error("Please select at least one source category.")
                    else:
                        # Create output directory
                        output_dir = Path(st.session_state.temp_dir) / "output"
                        output_dir.mkdir(exist_ok=True)

                        progress_bar = st.progress(0)
                        status_text = st.empty()

                        # Process each sample
                        sink_files = []
                        total_samples = len(pairs)

                        for idx, (base, info) in enumerate(pairs.items()):
                            sample_name = info['sample_name']
                            status_text.text(f"Processing {sample_name} ({idx+1}/{total_samples})...")

                            # Step 1: Extract V4 region from amplicons
                            with st.spinner(f"Extracting V4 region from {sample_name}..."):
                                trimmed_r1, trimmed_r2, trim_result = run_v4_extraction(
                                    info['R1'], info['R2'], output_dir, sample_name, selected_amplicon
                                )

                            if trimmed_r1 is None:
                                st.error(f"Error extracting V4 from {sample_name}:")
                                st.code(trim_result.stderr if trim_result else "Output files not created")
                                continue

                            # Step 2: Run DADA2 on trimmed files
                            spinner_msg = f"Running DADA2 on {sample_name}..."
                            if learn_error_rates:
                                spinner_msg += " (learning error rates - this will be slow)"
                            with st.spinner(spinner_msg):
                                result = run_dada2_processing(
                                    trimmed_r1, trimmed_r2,
                                    output_dir, sample_name, learn_error_rates
                                )

                            if result.returncode != 0:
                                st.error(f"Error processing {sample_name}:")
                                st.code(result.stderr)
                                continue

                            # Parse mapping stats from output
                            for line in result.stdout.split('\n'):
                                if line.startswith('MAPPING_STATS:'):
                                    stats = line.replace('MAPPING_STATS:', '').split(',')
                                    st.info(f"{sample_name}: Mapped {stats[0]}/{stats[1]} ASVs ({float(stats[2])/float(stats[3])*100:.1f}% reads)")

                            sink_file = output_dir / f"{sample_name}_abundance.tsv"
                            if sink_file.exists():
                                sink_files.append((sink_file, sample_name))

                            progress_bar.progress((idx + 1) / (total_samples + 2))

                        if sink_files:
                            # Create combined table
                            status_text.text("Creating combined feature table...")
                            combined_table, n_source_samples = create_combined_table(
                                sink_files, output_dir, selected_sources
                            )
                            st.info(f"Using {n_source_samples} source samples from {len(selected_sources)} categories")

                            progress_bar.progress((total_samples + 1) / (total_samples + 2))

                            # Create mapping file
                            sink_sample_names = [name for _, name in sink_files]
                            mapping_file = create_mapping_file(
                                sink_sample_names, output_dir, selected_sources
                            )

                            # Run SourceTracker2
                            status_text.text("Running SourceTracker2 Gibbs sampling...")
                            with st.spinner("Running SourceTracker2 (this may take a while)..."):
                                st2_output, st2_result = run_sourcetracker2(
                                    combined_table, mapping_file, output_dir
                                )

                            progress_bar.progress(1.0)

                            if st2_result.returncode == 0:
                                st.success("✅ Analysis complete!")
                                st.session_state.results_dir = st2_output
                                st.session_state.sink_samples = sink_sample_names
                                status_text.text("Analysis complete!")

                                # Show results inline
                                proportions, stds = load_results(st2_output)
                                if proportions is not None:
                                    st.divider()
                                    st.subheader("Source Contributions")
                                    for sample in proportions.index:
                                        st.plotly_chart(
                                            plot_source_contributions(proportions, sample),
                                            width='stretch'
                                        )
                            else:
                                st.error("SourceTracker2 failed:")
                                st.code(st2_result.stderr)
                        else:
                            st.error("No samples were successfully processed.")
            else:
                st.warning("Could not pair files. Please ensure R1 and R2 files follow naming convention.")
                st.write("Uploaded files:")
                for f in uploaded_files:
                    st.write(f"  - {f.name}")

    with tab2:
        st.header("Help & Information")

        st.subheader("About SourceTracker2")
        st.markdown("""
        SourceTracker2 is a Bayesian approach for predicting the source of microbial communities
        based on their composition. It estimates the proportion of sequences in a sink sample
        that originated from each source environment.
        """)

        st.subheader("File Format")
        st.markdown("""
        **Input files should be:**
        - Paired-end FASTQ files (R1 and R2)
        - Compressed with gzip (.fastq.gz)
        - Named with consistent pattern: `samplename_R1.fastq.gz` and `samplename_R2.fastq.gz`

        **Supported naming patterns:**
        - `sample_R1.fastq.gz` / `sample_R2.fastq.gz`
        - `sample_1.fastq.gz` / `sample_2.fastq.gz`
        """)

        st.subheader("Source Categories")
        st.markdown("""
        The database contains samples from the following source categories:
        """)

        source_info = get_source_counts(load_design_file())
        source_df = pd.DataFrame({
            'Category': list(source_info.keys()),
            'Sample Count': list(source_info.values())
        }).sort_values('Sample Count', ascending=False)
        st.dataframe(source_df, width='stretch', hide_index=True)

        st.subheader("Pipeline Steps")
        st.markdown("""
        1. **DADA2 Processing**: Filter, denoise, and merge paired-end reads
        2. **ASV Mapping**: Map sink ASVs to database ASVs by exact sequence match
        3. **SourceTracker2**: Run Gibbs sampling to estimate source proportions
        """)

if __name__ == "__main__":
    main()
