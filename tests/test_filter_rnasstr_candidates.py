import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "dataset"
    / "filter_rnasstr_candidates.py"
)


class FilterCandidatesTest(unittest.TestCase):
    def test_documented_filters_and_overlap_resolution(self):
        fields = [
            "id",
            "family",
            "sequence",
            "structure",
            "accession",
            "strand",
            "start",
            "end",
            "evalue",
        ]
        reference_rows = [
            {"id": "r1", "family": "RF00001", "sequence": "AACCCCGGUU", "structure": "((......))"},
            {"id": "r2", "family": "RF00001", "sequence": "AAGGCCCCUUUU", "structure": "((........))"},
            {"id": "r3", "family": "RF00001", "sequence": "AAGGGGCCCCUUUU", "structure": "((..........))"},
        ]
        candidate_rows = [
            {"id": "good", "family": "RF00001", "sequence": "AAGGCCCCUUUU", "structure": "((........))", "accession": "chr1", "strand": "+", "start": "1", "end": "12", "evalue": "1e-8"},
            {"id": "overlap", "family": "RF00001", "sequence": "GGAAAACCCCUU", "structure": "((........))", "accession": "chr1", "strand": "+", "start": "5", "end": "16", "evalue": "1e-4"},
            {"id": "opposite", "family": "RF00001", "sequence": "GGAAAACCCCUU", "structure": "((........))", "accession": "chr1", "strand": "-", "start": "5", "end": "16", "evalue": "1e-3"},
            {"id": "duplicate", "family": "RF00001", "sequence": "AAGGCCCCUUUU", "structure": "((........))", "accession": "chr2", "strand": "+", "start": "1", "end": "12", "evalue": "1e-6"},
            {"id": "long", "family": "RF00001", "sequence": "AA" + "C" * 16 + "UU", "structure": "((................))", "accession": "chr3", "strand": "+", "start": "1", "end": "20", "evalue": "1e-6"},
            {"id": "unpaired", "family": "RF00001", "sequence": "AAGGCCCCUUUU", "structure": "............", "accession": "chr4", "strand": "+", "start": "1", "end": "12", "evalue": "1e-6"},
            {"id": "ambiguous", "family": "RF00001", "sequence": "AAGGCCCCUUUN", "structure": "((........))", "accession": "chr5", "strand": "+", "start": "1", "end": "12", "evalue": "1e-6"},
            {"id": "unreferenced", "family": "RF99999", "sequence": "AACCGGUU", "structure": "((....))", "accession": "chr6", "strand": "+", "start": "1", "end": "8", "evalue": "1e-6"},
        ]

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.csv"
            candidates = root / "candidates.csv"
            output = root / "filtered.csv"
            rejections = root / "rejections.csv"
            summary = root / "summary.csv"

            with reference.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(reference_rows)
            with candidates.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(candidate_rows)

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--candidates",
                    str(candidates),
                    "--rfam-reference",
                    str(reference),
                    "--output",
                    str(output),
                    "--rejections",
                    str(rejections),
                    "--reference-summary",
                    str(summary),
                ],
                check=True,
            )

            with output.open(newline="") as handle:
                retained = list(csv.DictReader(handle))
            with rejections.open(newline="") as handle:
                rejected = list(csv.DictReader(handle))

            self.assertEqual(
                [row["id"] for row in retained],
                ["good", "opposite", "duplicate", "unreferenced"],
            )
            reasons = {row["id"]: row["reason"] for row in rejected}
            self.assertEqual(reasons["overlap"], "overlapping_hit")
            self.assertEqual(reasons["long"], "length_outlier")
            self.assertEqual(reasons["unpaired"], "low_pair_count")
            self.assertEqual(reasons["ambiguous"], "non_aucg_sequence")


if __name__ == "__main__":
    unittest.main()
