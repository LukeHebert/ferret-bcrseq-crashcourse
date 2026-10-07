"""Small end-to-end tests for the crashcourse plotting scripts."""

from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import plot_cluster_cdr3_diversity as cdr3_plot
import plot_lineage_v_usage as v_usage_plot
import plot_processing_summary as processing_plot


def write_fastq(path: Path, sequences: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        for index, sequence in enumerate(sequences, start=1):
            handle.write(f"@read{index}/1\n{sequence}\n+\n{'I' * len(sequence)}\n")


class PlotScriptTests(unittest.TestCase):
    def test_v_usage_deduplicates_lineages_and_normalizes_alleles(self) -> None:
        data = pd.DataFrame(
            {
                "ClusterID": [2, 1, 1, 3, 4],
                "v_call": ["IGHV2-1*02", "IGHV1-2*01,IGHV1-3*01", "IGHV1-2*01", "IGHV2-1*01", None],
            }
        )
        summary = v_usage_plot.summarise_lineage_v_usage(data)
        self.assertEqual(summary.to_dict("records"), [
            {"v_gene": "IGHV2-1", "lineage_count": 2},
            {"v_gene": "IGHV1-2", "lineage_count": 1},
        ])

    def test_cdr3_diversity_orders_numerically_and_limits_rows(self) -> None:
        data = pd.DataFrame(
            {
                "ClusterID": [10, 10, 2, 2, 1, 1],
                "cdr3_aa": ["AAA", "AAB", "BBB", "BBB", "CCC", "DDD"],
            }
        )
        summary = cdr3_plot.summarise_cluster_cdr3_diversity(data)
        self.assertEqual(summary["ClusterID"].tolist(), [1, 10, 2])
        self.assertEqual(summary["unique_cdr3_aa"].tolist(), [2, 2, 1])
        with tempfile.TemporaryDirectory() as temporary:
            png_path, tsv_path = cdr3_plot.write_outputs(summary, Path(temporary), top_n=2)
            self.assertTrue(png_path.is_file() and png_path.stat().st_size > 0)
            self.assertEqual(len(pd.read_csv(tsv_path, sep="\t")), 2)

    def test_processing_summary_counts_inputs_and_parses_filter_log(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            raw_r1 = directory / "raw_R1.fastq.gz"
            merged = directory / "merged.fastq.gz"
            write_fastq(raw_r1, ["ACGT", "TGCA", "GGGG"])
            write_fastq(merged, ["ACGT", "TGCA"])

            annotation = directory / "annotation.tsv"
            pd.DataFrame({"sequence": ["A", "B"]}).to_csv(annotation, sep="\t", index=False)
            filtered = directory / "filtered.tsv"
            pd.DataFrame({"sequence": ["A", "B"]}).to_csv(filtered, sep="\t", index=False)
            clustered = directory / "clustered.tsv"
            pd.DataFrame(
                {"sequence": ["A", "B"], "cdr3_aa": ["CAR", "CAR"], "ClusterID": [1, 2]}
            ).to_csv(clustered, sep="\t", index=False)
            log_directory = directory / "filter_collapse"
            log_directory.mkdir()
            (log_directory / "log_filtering_2026-01-01_00-00-00.txt").write_text(
                "Initial number of rows: 2\n"
                "After removing rows with NaN 'sequence': 2 rows left\n"
                "After removing rows where 'stop_codon' is 'T': 2 rows left\n"
                "Row count after removing duplicate sequence values: 2\n"
                "Row count after filtering sequences with count < 1: 2\n"
            )

            summary = processing_plot.make_summary(
                processing_plot.count_fastq_records(raw_r1),
                processing_plot.count_fastq_records(merged),
                len(processing_plot.read_nonempty_tsv(annotation)),
                processing_plot.parse_filter_stages(processing_plot.choose_filter_log(None, log_directory)),
                processing_plot.read_nonempty_tsv(filtered, processing_plot.FILTERED_REQUIRED),
                processing_plot.read_nonempty_tsv(clustered, processing_plot.CLUSTERED_REQUIRED),
            )
            counts = dict(zip(summary["category"], summary["count"]))
            self.assertEqual(counts["Raw paired reads"], 3)
            self.assertEqual(counts["Merged reads"], 2)
            self.assertEqual(counts["After removing duplicate sequence values"], 2)
            self.assertEqual(counts["Unique CDRH3 amino-acid sequences"], 1)
            self.assertEqual(counts["Unique ClusterID lineages"], 2)
            png_path, tsv_path = processing_plot.write_outputs(summary, directory / "plots")
            self.assertTrue(png_path.is_file() and png_path.stat().st_size > 0)
            self.assertTrue(tsv_path.is_file() and tsv_path.stat().st_size > 0)


if __name__ == "__main__":
    unittest.main()
