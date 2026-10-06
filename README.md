# Ferret BCR-seq crashcourse

This 90-minute, command-line workshop follows a small paired-end ferret BCR
dataset through quality control, read merging, IgBLAST annotation, sequence
collapse, and clonal clustering. The 5,000 read pairs are a teaching sample,
not a biologically representative dataset.

## Start a Codespace

1. On GitHub, select **Use this template** and create your own repository.
2. In that repository select **Code**, open the **Codespaces** tab, then select
   **Create codespace on main**. Choose the default two-core machine.
3. Wait for the setup notification. The container installs all tools and runs
   `make verify` automatically. In the terminal, run it again at any time:

   ```bash
   make verify
   ```

The environment contains FastQC, PEAR, cutadapt, Java, Python packages, and
IgBLAST 1.21.0 with the workshop's ferret reference overlay. You do not need
to install software on your own computer.

## What is in the demo data?

`data/ferret_demo/ferret_demo_R1.fastq.gz` and `_R2.fastq.gz` are a seeded,
uniform 5,000-pair selection from public SRA run `SRR33336518`. Its source
export had identical mate headers, so the teaching files append standard `/1`
and `/2` suffixes for PEAR compatibility. The pair count, seed, source export
command, and checksums are recorded in
`data/ferret_demo/MANIFEST.json`.

The original full FASTQs are intentionally not in this repository. To recreate
the sample after downloading `SRR33336518_1.fastq` and `_2.fastq` with
`fasterq-dump --split-files`, run:

```bash
python3 scripts/build_demo_data.py SRR33336518_1.fastq SRR33336518_2.fastq
```

## Run the workflow

Run each command from the repository root. Outputs are written beside the demo
inputs and are ignored by Git.

### 1. Quality-check, trim, and merge paired reads

```bash
python workflows/bcrseq_transcript/trim_merge.py \
  data/ferret_demo/ferret_demo_R1.fastq.gz \
  data/ferret_demo/ferret_demo_R2.fastq.gz \
  --threads 2
```

Inspect `data/ferret_demo/trim_merge/` for FastQC reports, PEAR logs, and
unassembled reads. The main output is `data/ferret_demo/ferret_demo.assembled.fastq`.

### 2. Annotate assembled reads with IgBLAST

```bash
python workflows/bcrseq_transcript/identify_genes.py \
  data/ferret_demo/ferret_demo.assembled.fastq ferret \
  --igblast-dir "$IGBLAST_HOME" --threads 2
```

This creates an AIRR-format table,
`data/ferret_demo/ferret_demo.assembled_IgBLAST.tsv`. Read the mapping log in
`data/ferret_demo/identify_genes/` to see the exact IgBLAST command.

### 3. Filter nonfunctional reads and collapse identical sequences

```bash
python workflows/bcrseq_transcript/filter_collapse.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST.tsv --threshold 1
```

The `_filtered.tsv` output retains functional sequences observed at least once;
`nt_seq_count` records how many identical nucleotide reads were collapsed. A
threshold of one is intentional here: this small, uniform teaching subset did
not retain repeated functional nucleotide reads at threshold two.

### 4. Cluster related heavy-chain sequences

```bash
python workflows/bcrseq_transcript/gupta_cluster.py \
  data/ferret_demo/ferret_demo.assembled_IgBLAST_filtered.tsv --auto_threshold
```

The final `_clustered.tsv` adds `ClusterID`. The `clustering/` directory
contains the selected threshold, a log, and a distance-to-nearest diagnostic
plot. Interpret cluster IDs only within this small demonstration sample.

## Instructor setup

Push this repository to GitHub, then mark it as a **template repository** in
its GitHub settings. After confirming a Codespace build succeeds, configure a
prebuild for `main` using `.devcontainer/devcontainer.json`; this makes student
creation faster but incurs Codespaces storage/compute usage. Students need
GitHub accounts with Codespaces access.

## Attribution

The workshop adapts the BCR transcript scripts under the repository's MIT
license. See `LICENSE` for the original copyright notice.
