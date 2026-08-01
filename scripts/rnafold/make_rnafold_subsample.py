#!/usr/bin/env python3
"""Create a deterministic family-stratified RNAfold benchmark subset."""

from __future__ import annotations

import argparse
import csv
import hashlib
import random
from collections import defaultdict
from pathlib import Path


def family_seed(global_seed: int, family: str) -> int:
    digest = hashlib.sha256(f"{global_seed}:{family}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--per-family", type=int, default=25)
    parser.add_argument("--max-length", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260724)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    all_counts = defaultdict(int)
    canonical_counts = defaultdict(int)
    eligible = defaultdict(list)
    with args.test.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"id", "sequence", "structure", "base_pairs", "len", "family"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"test CSV missing columns: {sorted(missing)}")
        for row in reader:
            family = row["family"]
            sequence = row["sequence"].upper().replace("T", "U")
            all_counts[family] += 1
            if set(sequence) <= set("ACGU"):
                canonical_counts[family] += 1
                if len(sequence) <= args.max_length:
                    row["sequence"] = sequence
                    eligible[family].append(row)

    selected = []
    for family in sorted(all_counts):
        candidates = eligible[family]
        rng = random.Random(family_seed(args.seed, family))
        if len(candidates) <= args.per_family:
            chosen = list(candidates)
        else:
            chosen = rng.sample(candidates, args.per_family)
        # Stable output ordering after sampling.
        chosen.sort(key=lambda row: row["id"])
        selected.extend(chosen)

    sample_csv = args.outdir / "rnafold_subsample.csv"
    sample_fasta = args.outdir / "rnafold_subsample.fasta"
    family_manifest = args.outdir / "rnafold_subsample_per_family.csv"
    run_manifest = args.outdir / "rnafold_subsample_manifest.txt"

    output_fields = [
        "rnafold_id", "id", "family", "sequence", "length",
        "ground_truth_dot_bracket", "ground_truth_base_pairs",
    ]
    with (
        sample_csv.open("w", newline="", encoding="utf-8") as csv_handle,
        sample_fasta.open("w", encoding="utf-8") as fasta_handle,
    ):
        writer = csv.DictWriter(csv_handle, fieldnames=output_fields)
        writer.writeheader()
        for index, row in enumerate(selected, 1):
            rnafold_id = f"rnafold_{index:06d}"
            writer.writerow({
                "rnafold_id": rnafold_id,
                "id": row["id"],
                "family": row["family"],
                "sequence": row["sequence"],
                "length": len(row["sequence"]),
                "ground_truth_dot_bracket": row["structure"],
                "ground_truth_base_pairs": row["base_pairs"],
            })
            fasta_handle.write(f">{rnafold_id}\n{row['sequence']}\n")

    manifest_fields = [
        "family", "test_sequences", "canonical_sequences",
        "eligible_sequences", "selected_sequences",
    ]
    with family_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields)
        writer.writeheader()
        for family in sorted(all_counts):
            writer.writerow({
                "family": family,
                "test_sequences": all_counts[family],
                "canonical_sequences": canonical_counts[family],
                "eligible_sequences": len(eligible[family]),
                "selected_sequences": min(len(eligible[family]), args.per_family),
            })

    lengths = [len(row["sequence"]) for row in selected]
    run_manifest.write_text(
        "\n".join([
            f"source_test_csv\t{args.test}",
            f"random_seed\t{args.seed}",
            f"target_per_family\t{args.per_family}",
            f"maximum_sequence_length\t{args.max_length}",
            "allowed_alphabet\tACGU",
            f"test_families\t{len(all_counts)}",
            f"represented_families\t{len({row['family'] for row in selected})}",
            f"selected_sequences\t{len(selected)}",
            f"minimum_selected_length\t{min(lengths)}",
            f"maximum_selected_length\t{max(lengths)}",
        ]) + "\n",
        encoding="utf-8",
    )
    for path in (sample_csv, sample_fasta, family_manifest, run_manifest):
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
