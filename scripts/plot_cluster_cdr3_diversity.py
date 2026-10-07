#!/usr/bin/env python3
"""Plot the number of unique CDRH3 amino-acid sequences in each lineage."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
from plotnine import aes, element_text, geom_col, ggplot, labs, theme, theme_minimal


REQUIRED_COLUMNS = {"ClusterID", "cdr3_aa"}


def summarise_cluster_cdr3_diversity(data: pd.DataFrame) -> pd.DataFrame:
    """Count distinct non-null CDRH3 amino-acid sequences per numeric ClusterID."""
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"Input TSV is missing required column(s): {', '.join(sorted(missing))}")

    lineages = data.dropna(subset=["ClusterID"]).copy()
    if lineages.empty:
        raise ValueError("Input TSV has no non-null ClusterID values")
    lineages["_cluster_sort"] = pd.to_numeric(lineages["ClusterID"], errors="coerce")
    if lineages["_cluster_sort"].isna().any():
        raise ValueError("ClusterID values must be numeric")

    summary = (
        lineages.groupby(["ClusterID", "_cluster_sort"], as_index=False)["cdr3_aa"]
        .nunique(dropna=True)
        .rename(columns={"cdr3_aa": "unique_cdr3_aa"})
        .sort_values(["unique_cdr3_aa", "_cluster_sort"], ascending=[False, True], kind="stable")
        .drop(columns="_cluster_sort")
        .reset_index(drop=True)
    )
    if summary.empty:
        raise ValueError("No lineage CDRH3 values were found")
    return summary


def make_plot(summary: pd.DataFrame):
    ordered = summary.copy()
    ordered["ClusterID"] = pd.Categorical(
        ordered["ClusterID"].astype(str), categories=ordered["ClusterID"].astype(str), ordered=True
    )
    return (
        ggplot(ordered, aes(x="ClusterID", y="unique_cdr3_aa"))
        + geom_col(fill="#4daf4a")
        + labs(
            title="CDRH3 amino-acid diversity by lineage",
            x="ClusterID",
            y="Unique CDRH3 amino-acid sequences",
        )
        + theme_minimal()
        + theme(axis_text_x=element_text(rotation=90, ha="right"))
    )


def write_outputs(summary: pd.DataFrame, output_dir: Path, top_n: int) -> tuple[Path, Path]:
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise ValueError(f"Could not create output directory {output_dir}: {error}") from error

    plotted = summary if top_n == 0 else summary.head(top_n).copy()
    tsv_path = output_dir / "cluster_cdr3_diversity.tsv"
    png_path = output_dir / "cluster_cdr3_diversity.png"
    plotted.to_csv(tsv_path, sep="\t", index=False)
    make_plot(plotted).save(
        str(png_path),
        width=max(8, 0.35 * len(plotted)),
        height=5,
        dpi=300,
        verbose=False,
    )
    return png_path, tsv_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("clustered_tsv", type=Path, help="Clustered TSV containing ClusterID and cdr3_aa")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for PNG and TSV outputs")
    parser.add_argument(
        "--top-n",
        type=int,
        default=30,
        help="Number of lineages to plot, ordered by diversity (default: 30; 0 plots all)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.top_n < 0:
        raise SystemExit("plot_cluster_cdr3_diversity: --top-n must be zero or positive")
    try:
        data = pd.read_csv(args.clustered_tsv, sep="\t")
        if data.empty:
            raise ValueError("Input TSV is empty")
        summary = summarise_cluster_cdr3_diversity(data)
        png_path, tsv_path = write_outputs(summary, args.output_dir, args.top_n)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        raise SystemExit(f"plot_cluster_cdr3_diversity: {error}")
    print(f"Wrote {png_path} and {tsv_path}")


if __name__ == "__main__":
    main()
