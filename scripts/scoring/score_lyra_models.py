#!/usr/bin/env python3
"""Normalize and score Lyra-TransPred outputs against RNASSTR truth.

F1 is recomputed from each raw file's TP/FP/FN counts. MCC uses the same
all-possible-unordered-nucleotide-pairs definition as the SincFold analysis.
Displayed dot-bracket strings are parsed separately and receive a serialization
status; serialization defects do not overwrite the model's original counts.
"""

from __future__ import annotations

import argparse
import ast
import csv
import math
from collections import defaultdict
from pathlib import Path


def parse_truth_pairs(value: str) -> set[tuple[int, int]]:
    pairs = set()
    for item in ast.literal_eval(value):
        i, j = int(item[0]), int(item[1])
        if i > j:
            i, j = j, i
        pairs.add((i, j))
    return pairs


def pairs_text(pairs: set[tuple[int, int]]) -> str:
    return "[" + ",".join(f"[{i},{j}]" for i, j in sorted(pairs)) + "]"


def parse_display_structure(
    structure: str | None, length: int
) -> tuple[str, str]:
    """Return enumerated displayed pairs and a serialization status."""
    if not structure:
        return "", "missing"
    if len(structure) != length:
        return "", "length_mismatch"
    stacks = {"(": [], "[": [], "{": [], "<": []}
    closer = {")": "(", "]": "[", "}": "{", ">": "<"}
    pairs: set[tuple[int, int]] = set()
    unmatched_closes = 0
    invalid_chars = 0
    for i, char in enumerate(structure):
        if char in stacks:
            stacks[char].append(i)
        elif char in closer:
            opener = closer[char]
            if stacks[opener]:
                pairs.add((stacks[opener].pop(), i))
            else:
                unmatched_closes += 1
        elif char != ".":
            invalid_chars += 1
    unmatched_opens = sum(len(stack) for stack in stacks.values())
    if invalid_chars:
        status = "invalid_characters"
    elif unmatched_opens or unmatched_closes:
        status = f"unbalanced_open_{unmatched_opens}_close_{unmatched_closes}"
    else:
        status = "valid"
    return pairs_text(pairs), status


def score_from_counts(
    tp: int, fp: int, fn: int, length: int
) -> tuple[int, float, float, float, float]:
    candidates = length * (length - 1) // 2
    tn = candidates - tp - fp - fn
    if min(tp, fp, fn, tn) < 0:
        raise ValueError("invalid confusion counts")
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0
    denominator = math.sqrt(
        (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    )
    mcc = (tp * tn - fp * fn) / denominator if denominator else 0.0
    return tn, precision, recall, f1, mcc


def load_raw(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, int]]:
    rows: dict[str, dict[str, str]] = {}
    stats = defaultdict(int)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {
            "id", "sequence", "family", "length", "true_structure",
            "predicted_structure", "tp", "fp", "fn",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for line_number, row in enumerate(reader, 2):
            stats["raw_rows"] += 1
            identifier = (row.get("id") or "").strip()
            if not identifier:
                stats["blank_id_rows"] += 1
                continue
            if identifier in rows:
                stats["duplicate_ids"] += 1
                continue
            rows[identifier] = row
    return rows, dict(stats)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--transpred", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    raw = {}
    load_stats = {}
    raw["transpred"], load_stats["transpred"] = load_raw(args.transpred)

    per_sequence = args.outdir / "lyra_per_sequence_scores.csv"
    global_summary = args.outdir / "lyra_global_summary.csv"
    family_summary = args.outdir / "lyra_per_family_summary.csv"
    validation = args.outdir / "lyra_validation_report.txt"

    base_fields = [
        "id", "family", "sequence", "length",
        "ground_truth_dot_bracket", "ground_truth_base_pairs",
    ]
    model_fields = [
        "status", "evaluation_scope", "evaluation_length", "prediction_dot_bracket",
        "prediction_base_pairs_from_display", "prediction_serialization_status",
        "tp", "fp", "fn", "tn", "precision", "recall", "f1", "mcc",
        "source_f1", "source_precision", "source_recall",
    ]
    fields = base_fields + [
        f"{model}_{field}"
        for model in ("transpred",)
        for field in model_fields
    ]
    totals = defaultdict(lambda: defaultdict(float))
    seen_truth = set()
    issues: list[str] = []

    with (
        args.truth.open(newline="", encoding="utf-8-sig") as truth_handle,
        per_sequence.open("w", newline="", encoding="utf-8") as output_handle,
    ):
        reader = csv.DictReader(truth_handle)
        writer = csv.DictWriter(output_handle, fieldnames=fields)
        writer.writeheader()
        for truth_line, truth in enumerate(reader, 2):
            identifier = truth["id"]
            seen_truth.add(identifier)
            sequence = truth["sequence"]
            length = len(sequence)
            truth_pairs = parse_truth_pairs(truth["base_pairs"])
            output = {
                "id": identifier,
                "family": truth["family"],
                "sequence": sequence,
                "length": length,
                "ground_truth_dot_bracket": truth["structure"],
                "ground_truth_base_pairs": pairs_text(truth_pairs),
            }
            for model in ("transpred",):
                prefix = f"{model}_"
                bucket = totals[(truth["family"], model)]
                bucket["n_truth"] += 1
                source = raw[model].get(identifier)
                if source is None:
                    output[prefix + "status"] = "missing"
                    bucket["n_missing"] += 1
                    continue
                try:
                    source_sequence = source["sequence"]
                    source_length = int(source["length"])
                    source_truth = source["true_structure"]
                    if (
                        source_sequence == sequence
                        and source_length == length
                        and source_truth == truth["structure"]
                    ):
                        evaluation_scope = "full_length"
                    elif (
                        source_length < length
                        and len(source_sequence) == source_length
                        and sequence.startswith(source_sequence)
                        and source_truth == truth["structure"][:source_length]
                    ):
                        evaluation_scope = "truncated_prefix"
                    else:
                        raise ValueError(
                            "source sequence/structure does not match full truth or a prefix"
                        )
                    tp, fp, fn = int(source["tp"]), int(source["fp"]), int(source["fn"])
                    tn, precision, recall, f1, mcc = score_from_counts(
                        tp, fp, fn, source_length
                    )
                    displayed_pairs, serialization_status = parse_display_structure(
                        source.get("predicted_structure"), source_length
                    )
                    output.update({
                        prefix + "status": (
                            "scored" if evaluation_scope == "full_length"
                            else "scored_truncated_prefix"
                        ),
                        prefix + "evaluation_scope": evaluation_scope,
                        prefix + "evaluation_length": source_length,
                        prefix + "prediction_dot_bracket":
                            source.get("predicted_structure") or "",
                        prefix + "prediction_base_pairs_from_display": displayed_pairs,
                        prefix + "prediction_serialization_status":
                            serialization_status,
                        prefix + "tp": tp,
                        prefix + "fp": fp,
                        prefix + "fn": fn,
                        prefix + "tn": tn,
                        prefix + "precision": precision,
                        prefix + "recall": recall,
                        prefix + "f1": f1,
                        prefix + "mcc": mcc,
                        prefix + "source_f1": source.get("f1", ""),
                        prefix + "source_precision": source.get("precision", ""),
                        prefix + "source_recall": source.get("recall", ""),
                    })
                    bucket["n_scored"] += 1
                    bucket["f1_sum"] += f1
                    bucket["mcc_sum"] += mcc
                    bucket["source_f1_abs_diff_sum"] += abs(
                        f1 - float(source.get("f1") or 0)
                    )
                    bucket[f"scope_{evaluation_scope}"] += 1
                    bucket[f"serialization_{serialization_status}"] += 1
                except Exception as exc:
                    output[prefix + "status"] = "invalid"
                    bucket["n_invalid"] += 1
                    issues.append(f"{model}\t{identifier}\t{exc}")
            writer.writerow(output)

    family_rows = []
    for (family, model), values in sorted(totals.items()):
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
    with family_summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(family_rows)

    global_rows = []
    for model in ("transpred",):
        model_buckets = [
            values for (family, name), values in totals.items() if name == model
        ]
        n_truth = sum(int(v["n_truth"]) for v in model_buckets)
        n_scored = sum(int(v["n_scored"]) for v in model_buckets)
        global_rows.append({
            "model": model,
            "n_truth": n_truth,
            "n_scored": n_scored,
            "n_missing": sum(int(v["n_missing"]) for v in model_buckets),
            "n_invalid": sum(int(v["n_invalid"]) for v in model_buckets),
            "coverage": n_scored / n_truth if n_truth else 0.0,
            "mean_f1": (
                sum(v["f1_sum"] for v in model_buckets) / n_scored
                if n_scored else ""
            ),
            "mean_mcc": (
                sum(v["mcc_sum"] for v in model_buckets) / n_scored
                if n_scored else ""
            ),
        })
    global_fields = summary_fields[1:]
    with global_summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=global_fields)
        writer.writeheader()
        writer.writerows(global_rows)

    report = [f"truth_rows\t{len(seen_truth)}"]
    for model in ("transpred",):
        model_buckets = [
            values for (family, name), values in totals.items() if name == model
        ]
        report.extend([
            f"{model}_raw_rows\t{load_stats[model].get('raw_rows', 0)}",
            f"{model}_loaded_ids\t{len(raw[model])}",
            f"{model}_blank_id_rows\t{load_stats[model].get('blank_id_rows', 0)}",
            f"{model}_duplicate_ids\t{load_stats[model].get('duplicate_ids', 0)}",
            f"{model}_unmatched_ids\t{len(set(raw[model]) - seen_truth)}",
            f"{model}_valid_serializations\t"
            f"{sum(int(v['serialization_valid']) for v in model_buckets)}",
            f"{model}_nonvalid_serializations\t"
            f"{sum(int(v['n_scored'] - v['serialization_valid']) for v in model_buckets)}",
            f"{model}_full_length_scores\t"
            f"{sum(int(v['scope_full_length']) for v in model_buckets)}",
            f"{model}_truncated_prefix_scores\t"
            f"{sum(int(v['scope_truncated_prefix']) for v in model_buckets)}",
            f"{model}_mean_abs_source_f1_difference\t"
            f"{sum(v['source_f1_abs_diff_sum'] for v in model_buckets) / max(1, sum(int(v['n_scored']) for v in model_buckets))}",
        ])
    report.append(f"invalid_rows\t{len(issues)}")
    if issues:
        report.extend(["", "model\tid\terror", *issues])
    validation.write_text("\n".join(report) + "\n", encoding="utf-8")

    for path in (per_sequence, global_summary, family_summary, validation):
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
