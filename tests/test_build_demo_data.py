"""Tests for deterministic, explicit demo-data selections."""

from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_demo_data.py"


def write_fastq(path: Path, mate_sequences: list[str]) -> None:
    with path.open("w") as handle:
        for index, sequence in enumerate(mate_sequences, start=1):
            handle.write(f"@SRR33336518.{index} source\n{sequence}\n+\n{'I' * len(sequence)}\n")


class BuildDemoDataTests(unittest.TestCase):
    def test_selection_file_rebuilds_source_order_reproducibly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            r1, r2 = directory / "R1.fastq", directory / "R2.fastq"
            write_fastq(r1, ["AAAA", "CCCC", "GGGG"])
            write_fastq(r2, ["TTTT", "GGGG", "CCCC"])
            selection = directory / "selected_spots.tsv"
            selection.write_text(
                "spot_id\tcategory\tlineage_key\n"
                "SRR33336518.3\tbackground\t\n"
                "SRR33336518.1\tenriched\tlineage-a\n"
            )

            outputs = []
            for name in ("first", "second"):
                out_dir = directory / name
                subprocess.run(
                    [
                        sys.executable, str(BUILDER), str(r1), str(r2),
                        "--out-dir", str(out_dir), "--pairs", "2",
                        "--selection-file", str(selection),
                    ],
                    check=True,
                    stdout=subprocess.DEVNULL,
                )
                outputs.append(out_dir)

            first_r1 = outputs[0] / "ferret_demo_R1.fastq.gz"
            second_r1 = outputs[1] / "ferret_demo_R1.fastq.gz"
            self.assertEqual(first_r1.read_bytes(), second_r1.read_bytes())
            with gzip.open(first_r1, "rt") as handle:
                self.assertEqual(handle.readline().strip(), "@SRR33336518.1/1 source")
                handle.readline(); handle.readline(); handle.readline()
                self.assertEqual(handle.readline().strip(), "@SRR33336518.3/1 source")
            manifest = json.loads((outputs[0] / "MANIFEST.json").read_text())
            self.assertEqual(manifest["sampling"]["categories"], {"background": 1, "enriched": 1})


if __name__ == "__main__":
    unittest.main()
