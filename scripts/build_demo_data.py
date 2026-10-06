#!/usr/bin/env python3
"""Build a deterministic paired FASTQ subset and provenance manifest.

The full SRA export is deliberately not tracked by the workshop repository.
Run this script with the two files made by ``fasterq-dump --split-files`` to
refresh the committed teaching subset.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
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
    args = p.parse_args()
    if args.pairs < 1:
        raise ValueError("--pairs must be positive")
    total = count_and_validate(args.r1, args.r2)
    if args.pairs > total:
        raise ValueError(f"Requested {args.pairs} pairs but source contains {total}")
    selected = set(random.Random(args.seed).sample(range(total), args.pairs))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output_r1 = args.out_dir / "ferret_demo_R1.fastq.gz"
    output_r2 = args.out_dir / "ferret_demo_R2.fastq.gz"
    with ExitStack() as stack:
        in_r1 = records(stack.enter_context(open_text(args.r1, "r")), str(args.r1))
        in_r2 = records(stack.enter_context(open_text(args.r2, "r")), str(args.r2))
        out_r1 = stack.enter_context(gzip.open(output_r1, "wt"))
        out_r2 = stack.enter_context(gzip.open(output_r2, "wt"))
        for index, (left, right) in enumerate(zip(in_r1, in_r2, strict=True)):
            if index in selected:
                out_r1.write("\n".join(normalized_record(left, 1)) + "\n")
                out_r2.write("\n".join(normalized_record(right, 2)) + "\n")
    manifest = {
        "accession": "SRR33336518",
        "source_export": "fasterq-dump --split-files SRR33336518",
        "source_pairs": total,
        "selected_pairs": args.pairs,
        "sampling": "uniform sample without replacement, retained in source order",
        "header_normalization": "SRA spot headers were suffixed with /1 and /2 for PEAR compatibility",
        "seed": args.seed,
        "files": {
            output_r1.name: {"sha256": digest(output_r1), "bytes": output_r1.stat().st_size},
            output_r2.name: {"sha256": digest(output_r2), "bytes": output_r2.stat().st_size},
        },
    }
    (args.out_dir / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
