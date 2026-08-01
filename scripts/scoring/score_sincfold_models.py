#!/usr/bin/env python3
"""Score published and retrained SincFold predictions against RNASSTR truth.

The script accepts the current SincFold CSV format (``id_clean`` plus
``base_pairs_predict``) and writes:

1. one combined per-sequence CSV containing truth and both model predictions;
2. one global summary CSV with coverage and macro-averaged F1/MCC per model;
3. one per-family summary CSV with coverage and macro-averaged F1/MCC.

Prediction base-pair indices are one-based in the inference files. Ground-truth
base-pair indices are zero-based. Output enumerated pairs are normalized to
zero-based indices so that all columns use one coordinate system.
"""

from __future__ import annotations

import argparse
import ast
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable


def sanitized_id(value: str) -> str:
    """Convert a ground-truth ID to the identifier used by SincFold output."""
    return value.replace("|", "_").replace("/", "_")


def parse_pairs(value: str | None, *, one_based: bool) -> set[tuple[int, int]]:
    """Parse an enumerated pair list and normalize it to zero-based tuples."""
    if value is None or not value.strip():
        return set()
    parsed = ast.literal_eval(value)
    if not isinstance(parsed, (list, tuple)):
        raise ValueError("base-pair value is not a list")

    pairs: set[tuple[int, int]] = set()
    for item in parsed:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            raise ValueError(f"invalid base-pair entry: {item!r}")
        i, j = int(item[0]), int(item[1])
        if one_based:
            i -= 1
            j -= 1
        if i == j:
            raise ValueError(f"self-pair is invalid: {(i, j)!r}")
        if i > j:
            i, j = j, i
        pairs.add((i, j))
    return pairs


def validate_pairs(pairs: Iterable[tuple[int, int]], length: int) -> None:
    for i, j in pairs:
        if i < 0 or j >= length:
            raise ValueError(f"pair {(i, j)} is outside sequence length {length}")


def pairs_text(pairs: set[tuple[int, int]]) -> str:
    return "[" + ",".join(f"[{i},{j}]" for i, j in sorted(pairs)) + "]"


def pairs_to_dot_bracket(pairs: set[tuple[int, int]], length: int) -> str:
    """Create dot-bracket notation; reject conflicting nucleotide partners."""
    chars = ["."] * length
    used: set[int] = set()
    for i, j in sorted(pairs):
        if i in used or j in used:
            raise ValueError("a nucleotide has multiple predicted partners")
        used.update((i, j))
        chars[i] = "("
        chars[j] = ")"
    return "".join(chars)


def metrics(
    truth: set[tuple[int, int]], prediction: set[tuple[int, int]], length: int
) -> tuple[int, int, int, int, float, float, float, float]:
    """Return exact-pair TP/FP/FN/TN, precision, recall, F1 and MCC.

    MCC treats each possible unordered nucleotide pair as a binary candidate;
    therefore there are length*(length-1)/2 candidates. This reproduces the
    convention in the existing SincFold scoring script and files.
    """
    tp = len(truth & prediction)
    fp = len(prediction - truth)
    fn = len(truth - prediction)
    candidates = length * (length - 1) // 2
    tn = candidates - tp - fp - fn
    if tn < 0:
        raise ValueError("negative TN count; pairs or sequence length are inconsistent")

    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0
    denominator = math.sqrt(
        (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    )
    mcc = (tp * tn - fp * fn) / denominator if denominator else 0.0
    return tp, fp, fn, tn, precision, recall, f1, mcc


def load_predictions(path: Path) -> tuple[dict[str, str], int]:
    """Load only identifiers and enumerated predictions from a SincFold CSV."""
    predictions: dict[str, str] = {}
    duplicate_count = 0
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"id_clean", "base_pairs_predict"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            identifier = (row["id_clean"] or "").strip()
            if not identifier:
                raise ValueError(f"{path}:{line_number}: blank id_clean")
            if identifier in predictions:
                duplicate_count += 1
                continue
            predictions[identifier] = row["base_pairs_predict"] or "[]"
    return predictions, duplicate_count


def add_summary(
    summary: dict[tuple[str, str], dict[str, float]],
    family: str,
    model: str,
    status: str,
    f1: float | None,
    mcc: float | None,
) -> None:
    bucket = summary[(family, model)]
    bucket["n_truth"] += 1
    if status == "scored":
        bucket["n_scored"] += 1
        bucket["f1_sum"] += float(f1)
        bucket["mcc_sum"] += float(mcc)
    elif status == "missing":
        bucket["n_missing"] += 1
    else:
        bucket["n_invalid"] += 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--published", required=True, type=Path)
    parser.add_argument("--retrained", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    model_paths = {
        "published": args.published,
        "retrained": args.retrained,
    }
    prediction_maps: dict[str, dict[str, str]] = {}
    duplicate_counts: dict[str, int] = {}
    for model, path in model_paths.items():
        prediction_maps[model], duplicate_counts[model] = load_predictions(path)

    per_sequence_path = args.outdir / "sincfold_per_sequence_scores.csv"
    global_path = args.outdir / "sincfold_global_summary.csv"
    family_path = args.outdir / "sincfold_per_family_summary.csv"
    validation_path = args.outdir / "sincfold_validation_report.txt"

    base_fields = [
        "id", "id_clean", "family", "sequence", "length",
        "ground_truth_dot_bracket", "ground_truth_base_pairs",
    ]
    metric_fields = [
        "status", "prediction_dot_bracket", "prediction_base_pairs",
        "tp", "fp", "fn", "tn", "precision", "recall", "f1", "mcc",
    ]
    output_fields = base_fields + [
        f"{model}_{field}"
        for model in ("published", "retrained")
        for field in metric_fields
    ]

    summary: dict[tuple[str, str], dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    seen_truth: set[str] = set()
    invalid_messages: list[str] = []

    with (
        args.truth.open(newline="", encoding="utf-8-sig") as truth_handle,
        per_sequence_path.open("w", newline="", encoding="utf-8") as out_handle,
    ):
        reader = csv.DictReader(truth_handle)
        required_truth = {"id", "sequence", "structure", "base_pairs", "family"}
        missing_truth = required_truth - set(reader.fieldnames or [])
        if missing_truth:
            raise ValueError(f"truth file missing columns {sorted(missing_truth)}")
        writer = csv.DictWriter(out_handle, fieldnames=output_fields)
        writer.writeheader()

        for line_number, row in enumerate(reader, start=2):
            identifier = row["id"]
            clean_id = sanitized_id(identifier)
            if clean_id in seen_truth:
                raise ValueError(f"duplicate truth identifier at line {line_number}: {clean_id}")
            seen_truth.add(clean_id)
            sequence = row["sequence"]
            length = len(sequence)
            if row.get("len") and int(row["len"]) != length:
                raise ValueError(f"truth length mismatch at line {line_number}: {identifier}")
            truth_pairs = parse_pairs(row["base_pairs"], one_based=False)
            validate_pairs(truth_pairs, length)

            output = {
                "id": identifier,
                "id_clean": clean_id,
                "family": row["family"],
                "sequence": sequence,
                "length": length,
                "ground_truth_dot_bracket": row["structure"],
                "ground_truth_base_pairs": pairs_text(truth_pairs),
            }

            for model, predictions in prediction_maps.items():
                raw_prediction = predictions.get(clean_id)
                prefix = f"{model}_"
                if raw_prediction is None:
                    output[prefix + "status"] = "missing"
                    add_summary(summary, row["family"], model, "missing", None, None)
                    continue
                try:
                    pred_pairs = parse_pairs(raw_prediction, one_based=True)
                    validate_pairs(pred_pairs, length)
                    pred_structure = pairs_to_dot_bracket(pred_pairs, length)
                    tp, fp, fn, tn, precision, recall, f1, mcc = metrics(
                        truth_pairs, pred_pairs, length
                    )
                    values = {
                        "status": "scored",
                        "prediction_dot_bracket": pred_structure,
                        "prediction_base_pairs": pairs_text(pred_pairs),
                        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                        "precision": precision, "recall": recall,
                        "f1": f1, "mcc": mcc,
                    }
                    output.update({prefix + key: value for key, value in values.items()})
                    add_summary(summary, row["family"], model, "scored", f1, mcc)
                except Exception as exc:  # retain the row and report malformed predictions
                    output[prefix + "status"] = "invalid"
                    output[prefix + "prediction_base_pairs"] = raw_prediction
                    add_summary(summary, row["family"], model, "invalid", None, None)
                    invalid_messages.append(f"{model}\t{clean_id}\t{exc}")
            writer.writerow(output)

    family_rows = []
    for (family, model), values in sorted(summary.items()):
        n_truth = int(values["n_truth"])
        n_scored = int(values["n_scored"])
        family_rows.append({
            "family": family,
            "model": model,
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

    global_rows = []
    for model in ("published", "retrained"):
        rows = [row for row in family_rows if row["model"] == model]
        n_truth = sum(int(row["n_truth"]) for row in rows)
        n_scored = sum(int(row["n_scored"]) for row in rows)
        # Reconstruct weighted sums from family means to preserve sequence-level macro means.
        f1_sum = sum(float(row["mean_f1"]) * int(row["n_scored"]) for row in rows if row["mean_f1"] != "")
        mcc_sum = sum(float(row["mean_mcc"]) * int(row["n_scored"]) for row in rows if row["mean_mcc"] != "")
        global_rows.append({
            "model": model,
            "n_truth": n_truth,
            "n_scored": n_scored,
            "n_missing": sum(int(row["n_missing"]) for row in rows),
            "n_invalid": sum(int(row["n_invalid"]) for row in rows),
            "coverage": n_scored / n_truth if n_truth else 0.0,
            "mean_f1": f1_sum / n_scored if n_scored else "",
            "mean_mcc": mcc_sum / n_scored if n_scored else "",
        })
    global_fields = [field for field in summary_fields if field != "family"]
    with global_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=global_fields)
        writer.writeheader()
        writer.writerows(global_rows)

    report_lines = [
        f"truth_rows\t{len(seen_truth)}",
        *(f"{model}_prediction_rows\t{len(prediction_maps[model])}" for model in prediction_maps),
        *(f"{model}_duplicate_ids\t{duplicate_counts[model]}" for model in prediction_maps),
        *(f"{model}_unmatched_prediction_ids\t{len(set(prediction_maps[model]) - seen_truth)}" for model in prediction_maps),
        f"invalid_predictions\t{len(invalid_messages)}",
    ]
    if invalid_messages:
        report_lines.extend(["", "model\tid_clean\terror", *invalid_messages])
    validation_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"Wrote {per_sequence_path}")
    print(f"Wrote {global_path}")
    print(f"Wrote {family_path}")
    print(f"Wrote {validation_path}")


if __name__ == "__main__":
    main()
