#!/usr/bin/env python3
"""Verify the workshop environment and the committed paired demo data."""

from __future__ import annotations

import gzip
import hashlib
import importlib
import json
import os
import shutil
import subprocess
import sys
import csv
from itertools import zip_longest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "ferret_demo"


def fastq_records(path: Path):
    with gzip.open(path, "rt") as handle:
        while True:
            lines = tuple(handle.readline().rstrip("\n") for _ in range(4))
            if not lines[0]:
                if any(lines[1:]):
                    raise AssertionError(f"{path}: truncated record")
                return
            assert lines[0].startswith("@") and lines[2].startswith("+"), f"{path}: malformed record"
            assert len(lines[1]) == len(lines[3]) == 300, f"{path}: expected 300 nt reads"
            yield lines


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pair_key(header: str, mate: str) -> str:
    token = header.split(maxsplit=1)[0]
    assert token.endswith(mate), f"Expected {mate} header suffix: {header}"
    return token[:-len(mate)]


def check_environment() -> None:
    for module in ("Bio", "pandas", "numpy", "scipy", "matplotlib"):
        importlib.import_module(module)
    for command in ("fastqc", "vsearch", "cutadapt"):
        assert shutil.which(command), f"Missing executable on PATH: {command}"
    home = Path(os.environ.get("IGBLAST_HOME", "/opt/ncbi-igblast-1.21.0"))
    igblastn = home / "bin" / "igblastn"
    assert igblastn.is_file(), f"Missing igblastn: {igblastn}"
    assert (home / "internal_data" / "ferret" / "ferret_V").is_file(), "Missing ferret V database"
    assert (home / "optional_file" / "ferret_gl.aux").is_file(), "Missing ferret auxiliary file"
    subprocess.run([str(igblastn), "-version"], check=True, stdout=subprocess.DEVNULL)


def check_data() -> None:
    manifest = json.loads((DATA / "MANIFEST.json").read_text())
    assert manifest["selected_pairs"] == 5000
    sampling = manifest.get("sampling")
    selection_spots: set[str] | None = None
    if isinstance(sampling, dict):
        selection_path = DATA / sampling["selection_file"]
        with selection_path.open(newline="") as handle:
            selection_rows = list(csv.DictReader(handle, delimiter="\t"))
        selection_spots = {row["spot_id"] for row in selection_rows}
        assert len(selection_spots) == len(selection_rows) == 5000, "Invalid selected-spot list"
        assert sha256(selection_path) == sampling["selection_file_sha256"], "Selection-list checksum mismatch"
        categories: dict[str, int] = {}
        for row in selection_rows:
            category = row.get("category", "selected") or "selected"
            categories[category] = categories.get(category, 0) + 1
        assert categories == sampling["categories"], "Selection category counts mismatch"
    curation = manifest.get("curation")
    if isinstance(curation, dict):
        lineages = curation["selected_lineages"]
        assert len(lineages) == 20, "Expected 20 selected full-data lineages"
        assert all(
            lineage["distinct_cdr3_aa"] >= 10 for lineage in lineages
        ), "Selected lineage has fewer than 10 distinct CDRH3 amino-acid members"
    pair_count = 0
    observed_spots: set[str] = set()
    for left, right in zip_longest(
        fastq_records(DATA / "ferret_demo_R1.fastq.gz"),
        fastq_records(DATA / "ferret_demo_R2.fastq.gz"),
    ):
        assert left is not None and right is not None, "R1/R2 have different record counts"
        assert pair_key(left[0], "/1") == pair_key(right[0], "/2"), (
            f"Header mismatch at pair {pair_count + 1}"
        )
        observed_spots.add(pair_key(left[0], "/1")[1:])
        pair_count += 1
    assert pair_count == 5000
    if selection_spots is not None:
        assert observed_spots == selection_spots, "FASTQ records differ from selected-spot list"
    for name, details in manifest["files"].items():
        assert sha256(DATA / name) == details["sha256"], f"Checksum mismatch: {name}"


def main() -> None:
    check_environment()
    check_data()
    print("Workshop environment and 5,000-pair demo data: OK")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, OSError, subprocess.CalledProcessError) as error:
        print(f"Workshop verification failed: {error}", file=sys.stderr)
        raise SystemExit(1)
