# Ferret BCR-seq crashcourse

This command-line crashcourse follows a 5,000-pair subset of μ-isotype-specific
ferret bulk BCR-seq paired-end reads through quality control, read merging,
IgBLAST annotation, sequence collapse, and clonal clustering. This teaching
sample is not intended for biological inference.

## Optional: discard an older workshop copy and start fresh

Use these steps only if you made an earlier personal copy of this workshop and
want to replace it with a fresh copy.

1. First, visit [your Codespaces](https://github.com/codespaces). Find the
   Codespace for the old workshop copy, click its **...** menu, and choose
   **Delete**. This removes its stored environment as well as avoiding future
   Codespaces storage use.
2. Open the old personal workshop repository on GitHub. Select **Settings**,
   remain on **General**, scroll to **Danger Zone**, and choose **Delete this
   repository**. Follow GitHub's confirmation prompts and type the repository
   name when asked.
3. Return to this template repository and follow the next section to create a
   new, up-to-date personal copy and Codespace.

## Start your own Codespace

These instructions assume that you are logged in to your GitHub account.

1. On this repository's GitHub page, click **Use this template** near the
   upper-right, above the file list. In its menu, click **Create a new
   repository**. Do not use the green **Code** button yet.
2. On the next page, choose your personal GitHub account as the owner, give the
   repository a name, and click **Create repository**. You are now viewing
   your own copy of the workshop.
3. In your copy, click the green **Code** button, select the **Codespaces**
   tab, then click **Create codespace on main**. Choose the default two-core
   machine.
4. Wait for the setup notification. The container installs all tools and runs
   `make verify` automatically. In the terminal, run it again at any time:

   ```bash
   make verify
   ```

5. When the command reports `Workshop environment and 5,000-pair demo data:
   OK`, continue with the workflow below.

## What is in the demo data?

`data/ferret_demo/ferret_demo_R1.fastq.gz` and
`data/ferret_demo/ferret_demo_R2.fastq.gz` contain 5,000 paired forward and
reverse reads from a μ-isotype-specific ferret bulk BCR-seq dataset. The pair
count and file checksums are recorded in `data/ferret_demo/MANIFEST.json`.

## Run the workflow

Run each command from the repository root. Outputs are written beside the demo
inputs and are ignored by Git.

### Run everything with a shell script

After you have worked through the individual steps below, inspect
`scripts/run_demo.sh` to see how a shell script can call each Python program
in sequence. `scripts/run_demo_literal.sh` contains the same workflow as a
plain, literal list of the commands shown below. Run either complete
demonstration with:

```bash
bash scripts/run_demo.sh
# or
bash scripts/run_demo_literal.sh
```

The script stops at the first error and writes the same outputs as the
individual commands, including the three plots.

### 1. Quality-check, trim, and merge paired reads

```bash
python workflows/bcrseq_transcript/trim_merge.py \
  data/ferret_demo/ferret_demo_R1.fastq.gz \
  data/ferret_demo/ferret_demo_R2.fastq.gz \
  --threads 2
```

Inspect `data/ferret_demo/trim_merge/` for FastQC reports and the VSEARCH merge
log. The main output is `data/ferret_demo/ferret_demo.assembled.fastq`.

To view the FastQC HTML reports in your browser, run this in a separate
terminal and leave it running:

```bash
python3 -m http.server 8000 --directory data/ferret_demo/trim_merge
```

Open the **Ports** tab, find port 8000, and click the globe/open-in-browser
icon. Stop the server with `Ctrl+C` before reusing port 8000 for another
directory.

### 2. Annotate assembled reads with IgBLAST

```bash
python workflows/bcrseq_transcript/identify_genes.py \
  data/ferret_demo/ferret_demo.assembled.fastq ferret \
  --igblast-dir "$IGBLAST_HOME" --threads 2
```

This creates an AIRR-format table,
`data/ferret_demo/ferret_demo.assembled_IgBLAST.tsv`. Read the mapping log in
`data/ferret_demo/identify_genes/` to see the exact IgBLAST command.

To browse that log through the same Codespaces browser route, run:

```bash
python3 -m http.server 8000 --directory data/ferret_demo/identify_genes
```

### 3. Filter nonfunctional reads and collapse identical sequences

```bash
python workflows/bcrseq_transcript/filter_collapse.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST.tsv --threshold 1
```

The `_filtered.tsv` output retains functional sequences observed at least once;
`nt_seq_count` records how many identical nucleotide reads were collapsed.
This teaching workflow uses a threshold of one so that every retained
functional nucleotide sequence is available for the clustering demonstration.

To view the filtering log in a browser, run:

```bash
python3 -m http.server 8000 --directory data/ferret_demo/filter_collapse
```

### 4. Cluster related heavy-chain sequences

```bash
python workflows/bcrseq_transcript/gupta_cluster.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered.tsv --auto_threshold
```

The final `_clustered.tsv` adds `ClusterID`. The `clustering/` directory
contains the selected threshold, a log, and a distance-to-nearest diagnostic
plot. Interpret cluster IDs only within this small demonstration sample.

To view the distance-to-nearest plot and clustering log in a browser, run:

```bash
python3 -m http.server 8000 --directory data/ferret_demo/clustering
```

## Visualize the results

These commands create PNG figures and tab-separated tables in
`data/ferret_demo/plots/`. Each bar chart is ordered from the largest value to
the smallest.

### V-gene usage among lineages

```bash
python scripts/plot_lineage_v_usage.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered_clustered.tsv \
  --output-dir data/ferret_demo/plots
```

This counts each `ClusterID` once and reports V genes without allele suffixes.

### CDRH3 amino-acid diversity within lineages

```bash
python scripts/plot_cluster_cdr3_diversity.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered_clustered.tsv \
  --output-dir data/ferret_demo/plots
```

The default figure shows the 30 lineages with the most unique CDRH3 amino-acid
sequences. Add `--top-n 0` to plot every lineage.

### Processing and diversity dashboard

```bash
python scripts/plot_processing_summary.py \
  --raw-r1 data/ferret_demo/ferret_demo_R1.fastq.gz \
  --merged-fastq data/ferret_demo/ferret_demo.assembled.fastq \
  --annotation-tsv data/ferret_demo/ferret_demo.assembled_IgBLAST.tsv \
  --filtered-tsv data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered.tsv \
  --clustered-tsv data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered_clustered.tsv \
  --filter-log-dir data/ferret_demo/filter_collapse \
  --output-dir data/ferret_demo/plots
```

The left panel shows raw, merged, annotated, and each logged filtering or
collapse count. The right panel summarizes final unique nucleotide sequences,
unique CDRH3 amino-acid sequences, and ClusterID lineages.

To view all three final PNG figures in a browser, run:

```bash
python3 -m http.server 8000 --directory data/ferret_demo/plots
```
