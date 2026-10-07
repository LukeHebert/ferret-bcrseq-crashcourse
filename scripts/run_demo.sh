#!/usr/bin/env bash
# Run the complete ferret BCR-seq crashcourse workflow and its three plots.
# Use: bash scripts/run_demo.sh

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python_bin="${PYTHON_BIN:-python3}"
if ! command -v "$python_bin" >/dev/null 2>&1; then
    echo "Python executable not found: $python_bin" >&2
    exit 1
fi
if [[ -z "${IGBLAST_HOME:-}" ]]; then
    echo "IGBLAST_HOME is not set. Open this project in its Codespace first." >&2
    exit 1
fi

data_dir="data/ferret_demo"
r1="$data_dir/ferret_demo_R1.fastq.gz"
r2="$data_dir/ferret_demo_R2.fastq.gz"
assembled="$data_dir/ferret_demo.assembled.fastq"
annotated="$data_dir/ferret_demo.assembled_IgBLAST.tsv"
filtered="$data_dir/ferret_demo.assembled_IgBLAST_filtered.tsv"
clustered="$data_dir/ferret_demo.assembled_IgBLAST_filtered_clustered.tsv"
plots="$data_dir/plots"

echo "[1/8] Checking the workshop environment and input data"
"$python_bin" scripts/verify_workshop.py

echo "[2/8] Quality-checking, trimming, and merging paired reads"
"$python_bin" workflows/bcrseq_transcript/trim_merge.py "$r1" "$r2" --threads 2

echo "[3/8] Annotating merged reads with IgBLAST"
"$python_bin" workflows/bcrseq_transcript/identify_genes.py \
    "$assembled" ferret --igblast-dir "$IGBLAST_HOME" --threads 2

echo "[4/8] Filtering and collapsing nucleotide sequences"
"$python_bin" workflows/bcrseq_transcript/filter_collapse.py "$annotated" --threshold 1

echo "[5/8] Clustering related heavy-chain sequences"
"$python_bin" workflows/bcrseq_transcript/gupta_cluster.py "$filtered" --auto_threshold

echo "[6/8] Plotting V-gene usage among lineages"
"$python_bin" scripts/plot_lineage_v_usage.py "$clustered" --output-dir "$plots"

echo "[7/8] Plotting CDRH3 diversity within lineages"
"$python_bin" scripts/plot_cluster_cdr3_diversity.py "$clustered" --output-dir "$plots"

echo "[8/8] Plotting the processing and diversity dashboard"
"$python_bin" scripts/plot_processing_summary.py \
    --raw-r1 "$r1" \
    --merged-fastq "$assembled" \
    --annotation-tsv "$annotated" \
    --filtered-tsv "$filtered" \
    --clustered-tsv "$clustered" \
    --filter-log-dir "$data_dir/filter_collapse" \
    --output-dir "$plots"

echo "Demo complete. Results and plots are in $data_dir."
