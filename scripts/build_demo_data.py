#!/usr/bin/env python3
"""Build a deterministic paired FASTQ subset and provenance manifest.

The full SRA export is deliberately not tracked by the workshop repository.
Run this script with the two files made by ``fasterq-dump --split-files`` to
refresh the committed teaching subset.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import random
from contextlib import ExitStack
from pathlib import Path
from typing import Iterator, TextIO


def open_text(path: Path, mode: str) -> TextIO:
    return gzip.open(path, mode + "t") if path.suffix == ".gz" else path.open(mode)


def records(handle: TextIO, label: str) -> Iterator[tuple[str, str, str, str]]:
    number = 0
    while True:
        lines = tuple(handle.readline().rstrip("\n") for _ in range(4))
        if not lines[0]:
            if any(lines[1:]):
                raise ValueError(f"{label}: truncated FASTQ record after record {number}")
            return
        if not all(lines) or not lines[0].startswith("@") or not lines[2].startswith("+"):
            raise ValueError(f"{label}: invalid FASTQ record {number + 1}")
        if len(lines[1]) != len(lines[3]):
            raise ValueError(f"{label}: sequence/quality length mismatch at record {number + 1}")
        number += 1
        yield lines  # type: ignore[misc]


def count_and_validate(r1: Path, r2: Path) -> int:
    total = 0
    with ExitStack() as stack:
        r1_records = records(stack.enter_context(open_text(r1, "r")), str(r1))
        r2_records = records(stack.enter_context(open_text(r2, "r")), str(r2))
        while True:
            left = next(r1_records, None)
            right = next(r2_records, None)
            if left is None or right is None:
                if left != right:
                    raise ValueError("R1 and R2 have different record counts")
                return total
            if left[0] != right[0]:
                raise ValueError(f"Paired headers differ at pair {total + 1}: {left[0]} != {right[0]}")
            total += 1


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def source_spot_id(header: str) -> str:
    """Return the SRA spot identifier from a raw FASTQ header."""
    return header[1:].split(maxsplit=1)[0].split("/", 1)[0]


def read_selection(path: Path) -> tuple[set[str], dict[str, int]]:
    """Read a committed spot-selection TSV and validate its required column."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "spot_id" not in rows[0]:
        raise ValueError(f"{path}: expected a nonempty TSV with a spot_id column")
    spots = [row["spot_id"] for row in rows]
    if any(not spot for spot in spots):
        raise ValueError(f"{path}: empty spot_id value")
    if len(spots) != len(set(spots)):
        raise ValueError(f"{path}: duplicate spot_id value")
    category_counts: dict[str, int] = {}
    for row in rows:
        category = row.get("category", "selected") or "selected"
        category_counts[category] = category_counts.get(category, 0) + 1
    return set(spots), category_counts


def deterministic_gzip_writer(stack: ExitStack, path: Path) -> TextIO:
    """Open a reproducible gzip text writer (zero timestamp and no filename)."""
    raw = stack.enter_context(path.open("wb"))
    compressed = stack.enter_context(gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename=""))
    return stack.enter_context(io.TextIOWrapper(compressed, newline="\n"))


def add_mate_suffix(line: str, mate: int) -> str:
    """Add a standard /1 or /2 suffix to an SRA header or plus line."""
    token, separator, remainder = line.partition(" ")
    return f"{token}/{mate}{separator}{remainder}"


def normalized_record(record: tuple[str, str, str, str], mate: int) -> tuple[str, str, str, str]:
    return (add_mate_suffix(record[0], mate), record[1],
            add_mate_suffix(record[2], mate), record[3])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("r1", type=Path)
    p.add_argument("r2", type=Path)
    p.add_argument("--out-dir", type=Path, default=Path("data/ferret_demo"))
    p.add_argument("--pairs", type=int, default=5000)
    p.add_argument("--seed", type=int, default=33336518)
    p.add_argument(
        "--selection-file", type=Path,
        help="TSV with a unique spot_id column; rebuild exactly those source pairs instead of random sampling.",
    )
    p.add_argument(
        "--curation-summary", type=Path,
        help="Optional JSON provenance summary to embed in MANIFEST.json for a curated selection.",
    )
    args = p.parse_args()
    if args.pairs < 1:
        raise ValueError("--pairs must be positive")
    total = count_and_validate(args.r1, args.r2)
    if args.pairs > total:
        raise ValueError(f"Requested {args.pairs} pairs but source contains {total}")
    category_counts: dict[str, int] | None = None
    if args.selection_file:
        selected_spots, category_counts = read_selection(args.selection_file)
        if len(selected_spots) != args.pairs:
            raise ValueError(
                f"{args.selection_file}: contains {len(selected_spots)} spots but --pairs is {args.pairs}"
            )
        selected: set[int] | set[str] = selected_spots
    else:
        selected = set(random.Random(args.seed).sample(range(total), args.pairs))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output_r1 = args.out_dir / "ferret_demo_R1.fastq.gz"
    output_r2 = args.out_dir / "ferret_demo_R2.fastq.gz"
    with ExitStack() as stack:
        in_r1 = records(stack.enter_context(open_text(args.r1, "r")), str(args.r1))
        in_r2 = records(stack.enter_context(open_text(args.r2, "r")), str(args.r2))
        out_r1 = deterministic_gzip_writer(stack, output_r1)
        out_r2 = deterministic_gzip_writer(stack, output_r2)
        for index, (left, right) in enumerate(zip(in_r1, in_r2, strict=True)):
            key = source_spot_id(left[0]) if args.selection_file else index
            if key in selected:
                out_r1.write("\n".join(normalized_record(left, 1)) + "\n")
                out_r2.write("\n".join(normalized_record(right, 2)) + "\n")
                if args.selection_file:
                    selected.remove(key)
    if args.selection_file and selected:
        preview = ", ".join(sorted(selected)[:3]) if args.selection_file else ""
        raise ValueError(f"Selection includes source spots not found in paired FASTQs: {preview}")
    if args.selection_file:
        sampling: str | dict[str, object] = {
            "method": "curated lineage-enriched selection, retained in source order",
            "selection_file": args.selection_file.name,
            "selection_file_sha256": digest(args.selection_file),
            "categories": category_counts,
            "seed": args.seed,
        }
    else:
        sampling = "uniform sample without replacement, retained in source order"
    manifest = {
        "accession": "SRR33336518",
        "source_export": "fasterq-dump --split-files SRR33336518",
        "source_pairs": total,
        "selected_pairs": args.pairs,
        "sampling": sampling,
        "header_normalization": "SRA spot headers were suffixed with /1 and /2 for standard paired-read compatibility",
        "seed": args.seed,
        "files": {
            output_r1.name: {"sha256": digest(output_r1), "bytes": output_r1.stat().st_size},
            output_r2.name: {"sha256": digest(output_r2), "bytes": output_r2.stat().st_size},
        },
    }
    if args.curation_summary:
        try:
            manifest["curation"] = json.loads(args.curation_summary.read_text())
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError(f"Could not read curation summary {args.curation_summary}: {error}") from error
    (args.out_dir / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
