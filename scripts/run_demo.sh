#!/usr/bin/env bash
# Run the complete workshop as the same literal commands shown in README.md.
# Use: bash scripts/run_demo.sh

set -e

python3 scripts/verify_workshop.py

python3 workflows/bcrseq_transcript/trim_merge.py \
  data/ferret_demo/ferret_demo_R1.fastq.gz \
  data/ferret_demo/ferret_demo_R2.fastq.gz \
  --threads 2

python3 workflows/bcrseq_transcript/identify_genes.py \
  data/ferret_demo/ferret_demo.assembled.fastq ferret \
  --igblast-dir "$IGBLAST_HOME" --threads 2

python3 workflows/bcrseq_transcript/filter_collapse.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST.tsv --threshold 1

python3 workflows/bcrseq_transcript/gupta_cluster.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered.tsv --auto_threshold

python3 scripts/plot_lineage_v_usage.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered_clustered.tsv \
  --output-dir data/ferret_demo/plots

python3 scripts/plot_cluster_cdr3_diversity.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered_clustered.tsv \
  --output-dir data/ferret_demo/plots

python3 scripts/plot_processing_summary.py \
  --raw-r1 data/ferret_demo/ferret_demo_R1.fastq.gz \
  --merged-fastq data/ferret_demo/ferret_demo.assembled.fastq \
  --annotation-tsv data/ferret_demo/ferret_demo.assembled_IgBLAST.tsv \
  --filtered-tsv data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered.tsv \
  --clustered-tsv data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered_clustered.tsv \
  --filter-log-dir data/ferret_demo/filter_collapse \
  --output-dir data/ferret_demo/plots
