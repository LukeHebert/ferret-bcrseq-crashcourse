#!/usr/bin/env python3
"""Plot V-gene usage after counting each ClusterID lineage once."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
from plotnine import aes, element_text, geom_col, ggplot, labs, theme, theme_minimal


REQUIRED_COLUMNS = {"ClusterID", "v_call"}


def v_gene(call: object) -> str | None:
    """Return the first assigned V gene without its allele suffix."""
    if pd.isna(call):
        return None
    first_call = str(call).split(",", 1)[0].strip()
    if not first_call:
        return None
    return first_call.split("*", 1)[0].strip() or None


def summarise_lineage_v_usage(data: pd.DataFrame) -> pd.DataFrame:
    """Count one normalized V-gene assignment per non-null ClusterID."""
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"Input TSV is missing required column(s): {', '.join(sorted(missing))}")

    lineages = data.dropna(subset=["ClusterID"]).drop_duplicates("ClusterID", keep="first").copy()
    if lineages.empty:
        raise ValueError("Input TSV has no non-null ClusterID values")
    lineages["v_gene"] = lineages["v_call"].map(v_gene)
    lineages = lineages.dropna(subset=["v_gene"])
    if lineages.empty:
        raise ValueError("No usable V-gene calls were found after lineage deduplication")

    return (
        lineages.groupby("v_gene", as_index=False)
        .size()
        .rename(columns={"size": "lineage_count"})
        .sort_values(["lineage_count", "v_gene"], ascending=[False, True], kind="stable")
        .reset_index(drop=True)
    )


def make_plot(summary: pd.DataFrame):
    ordered = summary.copy()
    ordered["v_gene"] = pd.Categorical(ordered["v_gene"], categories=ordered["v_gene"], ordered=True)
    return (
        ggplot(ordered, aes(x="v_gene", y="lineage_count"))
        + geom_col(fill="#377eb8")
        + labs(
            title="V-gene usage among unique lineages",
            x="V gene",
            y="Unique ClusterID lineages",
        )
        + theme_minimal()
        + theme(axis_text_x=element_text(rotation=45, ha="right"))
    )


def write_outputs(summary: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise ValueError(f"Could not create output directory {output_dir}: {error}") from error

    tsv_path = output_dir / "lineage_v_gene_usage.tsv"
    png_path = output_dir / "lineage_v_gene_usage.png"
    summary.to_csv(tsv_path, sep="\t", index=False)
    make_plot(summary).save(
        str(png_path),
        width=max(7, 0.45 * len(summary)),
        height=5,
        dpi=300,
        verbose=False,
    )
    return png_path, tsv_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("clustered_tsv", type=Path, help="Clustered TSV containing ClusterID and v_call")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for PNG and TSV outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        data = pd.read_csv(args.clustered_tsv, sep="\t")
        if data.empty:
            raise ValueError("Input TSV is empty")
        summary = summarise_lineage_v_usage(data)
        png_path, tsv_path = write_outputs(summary, args.output_dir)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        raise SystemExit(f"plot_lineage_v_usage: {error}")
    print(f"Wrote {png_path} and {tsv_path}")


if __name__ == "__main__":
    main()
