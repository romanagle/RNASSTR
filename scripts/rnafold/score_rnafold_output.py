#!/usr/bin/env python3
"""Score an RNAfold .fold output using the established RNASSTR metrics."""

from __future__ import annotations

import argparse
import ast
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


STRUCTURE_RE = re.compile(r"^([.()[\]{}<>]+)\s+\(\s*([-+0-9.eE]+)\s*\)")


def parse_pairs(value: str) -> set[tuple[int, int]]:
    pairs = set()
    for i, j in ast.literal_eval(value):
        i, j = int(i), int(j)
        if i > j:
            i, j = j, i
        pairs.add((i, j))
    return pairs


def dot_bracket_pairs(structure: str) -> set[tuple[int, int]]:
    stacks = {"(": [], "[": [], "{": [], "<": []}
    closer = {")": "(", "]": "[", "}": "{", ">": "<"}
    pairs = set()
    for index, char in enumerate(structure):
        if char in stacks:
            stacks[char].append(index)
        elif char in closer:
            opener = closer[char]
            if not stacks[opener]:
                raise ValueError("unmatched closing bracket")
            pairs.add((stacks[opener].pop(), index))
        elif char != ".":
            raise ValueError(f"invalid structure character {char!r}")
    if any(stacks.values()):
        raise ValueError("unmatched opening bracket")
    return pairs


def pairs_text(pairs: set[tuple[int, int]]) -> str:
    return "[" + ",".join(f"[{i},{j}]" for i, j in sorted(pairs)) + "]"


def metrics(truth, prediction, length):
    tp = len(truth & prediction)
    fp = len(prediction - truth)
    fn = len(truth - prediction)
    tn = length * (length - 1) // 2 - tp - fp - fn
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0
    denominator = math.sqrt(
        (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    )
    mcc = (tp * tn - fp * fn) / denominator if denominator else 0.0
    return tp, fp, fn, tn, precision, recall, f1, mcc


def read_rnafold(path: Path):
    predictions = {}
    current_id = None
    current_sequence = None
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                current_sequence = None
            elif current_id is not None and current_sequence is None:
                current_sequence = line.upper()
            elif current_id is not None:
                match = STRUCTURE_RE.match(line)
                if not match:
                    raise ValueError(
                        f"{path}:{line_number}: unexpected RNAfold result: {line}"
                    )
                if current_id in predictions:
                    raise ValueError(f"duplicate RNAfold ID: {current_id}")
                predictions[current_id] = {
                    "sequence": current_sequence,
                    "structure": match.group(1),
                    "mfe_kcal_mol": float(match.group(2)),
                }
                current_id = None
                current_sequence = None
    return predictions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, type=Path)
    parser.add_argument("--rnafold-output", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    predictions = read_rnafold(args.rnafold_output)

    per_sequence_path = args.outdir / "rnafold_per_sequence_scores.csv"
    global_path = args.outdir / "rnafold_global_summary.csv"
    family_path = args.outdir / "rnafold_per_family_summary.csv"
    validation_path = args.outdir / "rnafold_validation_report.txt"
    fields = [
        "rnafold_id", "id", "family", "sequence", "length",
        "ground_truth_dot_bracket", "ground_truth_base_pairs",
        "status", "prediction_dot_bracket", "prediction_base_pairs",
        "mfe_kcal_mol", "tp", "fp", "fn", "tn",
        "precision", "recall", "f1", "mcc",
    ]
    totals = defaultdict(lambda: defaultdict(float))
    seen = set()
    invalid = []
    with (
        args.sample.open(newline="", encoding="utf-8") as sample_handle,
        per_sequence_path.open("w", newline="", encoding="utf-8") as out_handle,
    ):
        reader = csv.DictReader(sample_handle)
        writer = csv.DictWriter(out_handle, fieldnames=fields)
        writer.writeheader()
        for row in reader:
            rnafold_id = row["rnafold_id"]
            seen.add(rnafold_id)
            family = row["family"]
            bucket = totals[family]
            bucket["n_truth"] += 1
            output = dict(row)
            prediction = predictions.get(rnafold_id)
            if prediction is None:
                output["status"] = "missing"
                bucket["n_missing"] += 1
            else:
                try:
                    if prediction["sequence"] != row["sequence"]:
                        raise ValueError("sequence mismatch")
                    if len(prediction["structure"]) != int(row["length"]):
                        raise ValueError("structure length mismatch")
                    truth_pairs = parse_pairs(row["ground_truth_base_pairs"])
                    pred_pairs = dot_bracket_pairs(prediction["structure"])
                    values = metrics(truth_pairs, pred_pairs, int(row["length"]))
                    output.update({
                        "status": "scored",
                        "prediction_dot_bracket": prediction["structure"],
                        "prediction_base_pairs": pairs_text(pred_pairs),
                        "mfe_kcal_mol": prediction["mfe_kcal_mol"],
                    })
                    for name, value in zip(
                        ("tp", "fp", "fn", "tn", "precision", "recall", "f1", "mcc"),
                        values,
                    ):
                        output[name] = value
                    bucket["n_scored"] += 1
                    bucket["f1_sum"] += values[-2]
                    bucket["mcc_sum"] += values[-1]
                except Exception as exc:
                    output["status"] = "invalid"
                    bucket["n_invalid"] += 1
                    invalid.append(f"{rnafold_id}\t{exc}")
            writer.writerow(output)

    family_rows = []
    for family, values in sorted(totals.items()):
        n_truth, n_scored = int(values["n_truth"]), int(values["n_scored"])
        family_rows.append({
            "family": family,
            "model": "rnafold",
            "n_truth": n_truth,
            "n_scored": n_scored,
            "n_missing": int(values["n_missing"]),
            "n_invalid": int(values["n_invalid"]),
            "coverage": n_scored / n_truth if n_truth else 0.0,
            "mean_f1": values["f1_sum"] / n_scored if n_scored else "",
            "mean_mcc": values["mcc_sum"] / n_scored if n_scored else "",
        })
    summary_fields = [
        "family", "model", "n_truth", "n_scored", "n_missing", "n_invalid",
        "coverage", "mean_f1", "mean_mcc",
    ]
    with family_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(family_rows)

    n_truth = sum(row["n_truth"] for row in family_rows)
    n_scored = sum(row["n_scored"] for row in family_rows)
    global_row = {
        "model": "rnafold",
        "n_truth": n_truth,
        "n_scored": n_scored,
        "n_missing": sum(row["n_missing"] for row in family_rows),
        "n_invalid": sum(row["n_invalid"] for row in family_rows),
        "coverage": n_scored / n_truth if n_truth else 0.0,
        "mean_f1": sum(
            float(row["mean_f1"]) * row["n_scored"]
            for row in family_rows if row["mean_f1"] != ""
        ) / n_scored if n_scored else "",
        "mean_mcc": sum(
            float(row["mean_mcc"]) * row["n_scored"]
            for row in family_rows if row["mean_mcc"] != ""
        ) / n_scored if n_scored else "",
    }
    with global_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields[1:])
        writer.writeheader()
        writer.writerow(global_row)

    unmatched = set(predictions) - seen
    report = [
        f"sample_rows\t{len(seen)}",
        f"rnafold_output_rows\t{len(predictions)}",
        f"unmatched_output_ids\t{len(unmatched)}",
        f"invalid_rows\t{len(invalid)}",
    ]
    if invalid:
        report.extend(["", "rnafold_id\terror", *invalid])
    validation_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    for path in (per_sequence_path, global_path, family_path, validation_path):
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
