#!/usr/bin/env python3
"""Create a two-panel summary of BCR-seq processing and final repertoire diversity."""

from __future__ import annotations

import argparse
import gzip
from io import BytesIO
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from plotnine import aes, coord_flip, geom_col, ggplot, labs, theme, theme_minimal


FILTER_PATTERNS = (
    re.compile(r"^After (?P<label>.+): (?P<count>\d+) rows left$"),
    re.compile(r"^Row count after (?P<label>.+): (?P<count>\d+)$"),
)
FILTERED_REQUIRED = {"sequence"}
CLUSTERED_REQUIRED = {"ClusterID", "cdr3_aa"}


def count_fastq_records(path: Path) -> int:
    """Validate a FASTQ/FASTQ.GZ file and return its record count."""
    opener = gzip.open if path.suffix == ".gz" else open
    count = 0
    try:
        with opener(path, "rt") as handle:
            while True:
                lines = [handle.readline().rstrip("\n") for _ in range(4)]
                if not lines[0]:
                    if any(lines[1:]):
                        raise ValueError("truncated record")
                    break
                if not lines[0].startswith("@") or not lines[2].startswith("+"):
                    raise ValueError("malformed record")
                if len(lines[1]) != len(lines[3]):
                    raise ValueError("sequence and quality lengths differ")
                count += 1
    except OSError as error:
        raise ValueError(f"Could not read FASTQ {path}: {error}") from error
    if count == 0:
        raise ValueError(f"FASTQ {path} contains no records")
    return count


def read_nonempty_tsv(path: Path, required_columns: set[str] | None = None) -> pd.DataFrame:
    try:
        data = pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.ParserError) as error:
        raise ValueError(f"Could not read TSV {path}: {error}") from error
    if data.empty:
        raise ValueError(f"TSV {path} contains no rows")
    if required_columns:
        missing = required_columns.difference(data.columns)
        if missing:
            raise ValueError(f"TSV {path} is missing required column(s): {', '.join(sorted(missing))}")
    return data


def choose_filter_log(log_path: Path | None, log_dir: Path | None) -> Path:
    if log_path is not None:
        if not log_path.is_file():
            raise ValueError(f"Filter log does not exist: {log_path}")
        return log_path
    assert log_dir is not None
    if not log_dir.is_dir():
        raise ValueError(f"Filter-log directory does not exist: {log_dir}")
    logs = sorted(log_dir.glob("log_filtering_*.txt"), key=lambda item: item.stat().st_mtime, reverse=True)
    if not logs:
        raise ValueError(f"No log_filtering_*.txt files found in {log_dir}")
    return logs[0]


def parse_filter_stages(log_path: Path) -> list[tuple[str, int]]:
    """Return each ordered filtering/collapse count from a filter_collapse log."""
    stages: list[tuple[str, int]] = []
    try:
        lines = log_path.read_text().splitlines()
    except OSError as error:
        raise ValueError(f"Could not read filter log {log_path}: {error}") from error
    for line in lines:
        for pattern in FILTER_PATTERNS:
            match = pattern.match(line)
            if match:
                stages.append((match.group("label"), int(match.group("count"))))
                break
    if not stages:
        raise ValueError(f"Filter log has no recognized row-count stages: {log_path}")
    return stages


def make_summary(
    raw_pairs: int,
    merged_reads: int,
    annotation_rows: int,
    filter_stages: list[tuple[str, int]],
    filtered: pd.DataFrame,
    clustered: pd.DataFrame,
) -> pd.DataFrame:
    process_rows = [
        {"panel": "Processing", "category": "Raw paired reads", "count": raw_pairs},
        {"panel": "Processing", "category": "Merged reads", "count": merged_reads},
        {"panel": "Processing", "category": "IgBLAST annotations", "count": annotation_rows},
    ]
    process_rows.extend(
        {"panel": "Processing", "category": f"After {label}", "count": count}
        for label, count in filter_stages
    )
    diversity_rows = [
        {
            "panel": "Final diversity",
            "category": "Filtered unique nucleotide sequences",
            "count": int(filtered["sequence"].nunique(dropna=True)),
        },
        {
            "panel": "Final diversity",
            "category": "Unique CDRH3 amino-acid sequences",
            "count": int(clustered["cdr3_aa"].nunique(dropna=True)),
        },
        {
            "panel": "Final diversity",
            "category": "Unique ClusterID lineages",
            "count": int(clustered["ClusterID"].nunique(dropna=True)),
        },
    ]
    return pd.DataFrame(process_rows + diversity_rows)


def horizontal_bar_plot(data: pd.DataFrame, title: str):
    ordered = data.copy()
    categories = list(reversed(ordered["category"].tolist()))
    ordered["category"] = pd.Categorical(ordered["category"], categories=categories, ordered=True)
    return (
        ggplot(ordered, aes(x="category", y="count"))
        + geom_col(fill="#984ea3")
        + coord_flip()
        + labs(title=title, x="", y="Count")
        + theme_minimal()
        + theme(figure_size=(8, max(3.5, 0.45 * len(ordered))))
    )


def write_outputs(summary: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise ValueError(f"Could not create output directory {output_dir}: {error}") from error
    tsv_path = output_dir / "processing_summary.tsv"
    png_path = output_dir / "processing_summary.png"
    summary.to_csv(tsv_path, sep="\t", index=False)
    processing = horizontal_bar_plot(summary[summary["panel"] == "Processing"], "Read processing")
    diversity = horizontal_bar_plot(summary[summary["panel"] == "Final diversity"], "Final repertoire diversity")
    try:
        # Plotnine 0.15+ supports native plot composition.
        dashboard = (processing | diversity) + theme(figure_size=(18, 6))
        dashboard.save(str(png_path), dpi=300)
    except TypeError:
        # Each panel remains a Plotnine plot. Matplotlib only arranges the
        # rendered panels for older local Plotnine installations.
        images = []
        for panel in (processing, diversity):
            figure = panel.draw()
            buffer = BytesIO()
            figure.savefig(buffer, format="png", dpi=300, bbox_inches="tight")
            plt.close(figure)
            buffer.seek(0)
            images.append(plt.imread(buffer))
        figure, axes = plt.subplots(1, 2, figsize=(18, 6))
        for axis, image in zip(axes, images, strict=True):
            axis.imshow(image)
            axis.axis("off")
        figure.savefig(png_path, dpi=300, bbox_inches="tight")
        plt.close(figure)
    return png_path, tsv_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-r1", required=True, type=Path, help="Raw R1 FASTQ or FASTQ.GZ")
    parser.add_argument("--merged-fastq", required=True, type=Path, help="Merged assembled FASTQ")
    parser.add_argument("--annotation-tsv", required=True, type=Path, help="IgBLAST annotation TSV")
    parser.add_argument("--filtered-tsv", required=True, type=Path, help="Filtered/collapsed TSV")
    parser.add_argument("--clustered-tsv", required=True, type=Path, help="Final clustered TSV")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for PNG and TSV outputs")
    logs = parser.add_mutually_exclusive_group(required=True)
    logs.add_argument("--filter-log", type=Path, help="Specific filter_collapse log file")
    logs.add_argument("--filter-log-dir", type=Path, help="Directory containing filter_collapse logs; newest is used")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        annotation = read_nonempty_tsv(args.annotation_tsv)
        filtered = read_nonempty_tsv(args.filtered_tsv, FILTERED_REQUIRED)
        clustered = read_nonempty_tsv(args.clustered_tsv, CLUSTERED_REQUIRED)
        filter_log = choose_filter_log(args.filter_log, args.filter_log_dir)
        summary = make_summary(
            count_fastq_records(args.raw_r1),
            count_fastq_records(args.merged_fastq),
            len(annotation),
            parse_filter_stages(filter_log),
            filtered,
            clustered,
        )
        png_path, tsv_path = write_outputs(summary, args.output_dir)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        raise SystemExit(f"plot_processing_summary: {error}")
    print(f"Wrote {png_path} and {tsv_path}")


if __name__ == "__main__":
    main()
