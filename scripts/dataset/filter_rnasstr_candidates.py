#!/usr/bin/env python3
"""Apply the documented RNASSTR family-level quality-control criteria."""

from __future__ import annotations

import argparse
import ast
import csv
import math
import sqlite3
import statistics
import tempfile
from collections import defaultdict
from pathlib import Path


CANONICAL_PAIRS = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),
}
VALID_BASES = set("AUCG")
OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}", "<": ">"}
OPEN_TO_CLOSE.update({chr(code): chr(code).lower() for code in range(65, 91)})
CLOSE_TO_OPEN = {close: open_ for open_, close in OPEN_TO_CLOSE.items()}


def family_value(row):
    for field in ("family", "rfam_id", "family_key"):
        value = (row.get(field) or "").strip()
        if value:
            return value
    raise ValueError("record lacks family/rfam_id/family_key")


def parse_pairs(row):
    value = (row.get("base_pairs") or "").strip()
    if value:
        parsed = ast.literal_eval(value)
        pairs = set()
        for item in parsed:
            if not isinstance(item, (list, tuple)) or len(item) != 2:
                raise ValueError(f"invalid base-pair entry: {item!r}")
            i, j = int(item[0]), int(item[1])
            if i > j:
                i, j = j, i
            if i == j:
                raise ValueError(f"self pair: {(i, j)!r}")
            pairs.add((i, j))
        return pairs

    structure = (row.get("structure") or "").strip()
    if not structure:
        raise ValueError("record lacks structure and base_pairs")
    stacks = defaultdict(list)
    pairs = set()
    for index, character in enumerate(structure):
        if character in OPEN_TO_CLOSE:
            stacks[character].append(index)
        elif character in CLOSE_TO_OPEN:
            opener = CLOSE_TO_OPEN[character]
            if stacks[opener]:
                pairs.add((stacks[opener].pop(), index))
    return pairs


def record_metrics(row):
    sequence = (row.get("sequence") or "").upper()
    if not sequence:
        raise ValueError("empty sequence")
    invalid_bases = sorted(set(sequence) - VALID_BASES)
    if invalid_bases:
        raise ValueError(f"sequence contains non-AUCG characters: {''.join(invalid_bases)}")
    pairs = parse_pairs(row)
    for i, j in pairs:
        if i < 0 or j >= len(sequence):
            raise ValueError(f"base pair {(i, j)} exceeds sequence length")
    canonical = sum(
        (sequence[i], sequence[j]) in CANONICAL_PAIRS for i, j in pairs
    )
    return sequence, len(sequence), len(pairs), canonical


def summary(values):
    mean = statistics.fmean(values)
    standard_deviation = statistics.pstdev(values) if len(values) > 1 else 0.0
    return mean, standard_deviation


def load_reference_statistics(path):
    measurements = defaultdict(lambda: {"length": [], "pairs": [], "canonical": []})
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            family = family_value(row)
            _, length, pair_count, canonical_count = record_metrics(row)
            measurements[family]["length"].append(length)
            measurements[family]["pairs"].append(pair_count)
            measurements[family]["canonical"].append(canonical_count)

    statistics_by_family = {}
    for family, values in measurements.items():
        length_mean, length_sd = summary(values["length"])
        pair_mean, pair_sd = summary(values["pairs"])
        canonical_mean, canonical_sd = summary(values["canonical"])
        statistics_by_family[family] = {
            "reference_n": len(values["length"]),
            "length_mean": length_mean,
            "length_sd": length_sd,
            "pair_mean": pair_mean,
            "pair_sd": pair_sd,
            "canonical_mean": canonical_mean,
            "canonical_sd": canonical_sd,
        }
    return statistics_by_family


def parse_optional_float(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.inf
    return result if math.isfinite(result) else math.inf


def parse_optional_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def interval_fields(row):
    accession = (row.get("accession") or "").strip()
    strand = (row.get("strand") or row.get("hit_strand") or "").strip()
    start = parse_optional_int(row.get("start"))
    end = parse_optional_int(row.get("end"))
    if not accession or start is None or end is None:
        return None
    return accession, strand, min(start, end), max(start, end)


def initialize_database(connection):
    connection.executescript(
        """
        CREATE TABLE intervals (
            row_number INTEGER PRIMARY KEY,
            accession TEXT NOT NULL,
            strand TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            evalue REAL NOT NULL,
            record_id TEXT NOT NULL
        );
        CREATE INDEX intervals_accession_start
            ON intervals(accession, strand, start, end);
        """
    )


def choose_interval_winner(component):
    return min(component, key=lambda row: row[5])


def resolve_overlaps(connection):
    rejected = {}
    cursor = connection.execute(
        """
        SELECT row_number, accession, strand, start, end, evalue, record_id
        FROM intervals
        ORDER BY accession, strand, start, end
        """
    )
    component = []
    component_key = None
    component_end = None

    def finish(rows):
        if len(rows) < 2:
            return
        winner = choose_interval_winner(rows)
        for row in rows:
            if row[0] != winner[0]:
                rejected[row[0]] = winner[6]

    for row in cursor:
        accession, strand, start, end = row[1], row[2], row[3], row[4]
        key = (accession, strand)
        if (
            component
            and key == component_key
            and start <= component_end
        ):
            component.append(row)
            component_end = max(component_end, end)
        else:
            finish(component)
            component = [row]
            component_key = key
            component_end = end
    finish(component)
    return rejected


def write_reference_summary(path, statistics_by_family, standard_deviations):
    fields = [
        "family",
        "reference_n",
        "length_mean",
        "length_sd",
        "length_min",
        "length_max",
        "pair_mean",
        "pair_sd",
        "pair_min",
        "canonical_mean",
        "canonical_sd",
        "canonical_min",
    ]
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for family, values in sorted(statistics_by_family.items()):
            writer.writerow({
                "family": family,
                **values,
                "length_min": max(0, values["length_mean"] - standard_deviations * values["length_sd"]),
                "length_max": values["length_mean"] + standard_deviations * values["length_sd"],
                "pair_min": max(0, values["pair_mean"] - standard_deviations * values["pair_sd"]),
                "canonical_min": max(0, values["canonical_mean"] - standard_deviations * values["canonical_sd"]),
            })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--rfam-reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--rejections", required=True, type=Path)
    parser.add_argument("--reference-summary", required=True, type=Path)
    parser.add_argument("--standard-deviations", type=float, default=2.0)
    parser.add_argument("--skip-overlap-resolution", action="store_true")
    args = parser.parse_args()

    reference = load_reference_statistics(args.rfam_reference)
    write_reference_summary(
        args.reference_summary, reference, args.standard_deviations
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.rejections.parent.mkdir(parents=True, exist_ok=True)
    args.reference_summary.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="rnasstr_qc_") as temporary_directory:
        temporary_path = Path(temporary_directory)
        database = sqlite3.connect(temporary_path / "qc.sqlite")
        initialize_database(database)

        with (
            args.candidates.open(newline="", encoding="utf-8-sig") as input_handle,
            (temporary_path / "accepted_pre_overlap.csv").open(
                "w", newline="", encoding="utf-8"
            ) as preliminary_handle,
            args.rejections.open("w", newline="", encoding="utf-8") as reject_handle,
        ):
            reader = csv.DictReader(input_handle)
            if not reader.fieldnames:
                raise ValueError("candidate CSV has no header")
            output_fields = reader.fieldnames + [
                "qc_sequence_length",
                "qc_pair_count",
                "qc_canonical_pair_count",
            ]
            preliminary_writer = csv.DictWriter(
                preliminary_handle, fieldnames=["qc_row_number"] + output_fields
            )
            preliminary_writer.writeheader()
            rejection_writer = csv.DictWriter(
                reject_handle,
                fieldnames=["id", "family", "reason", "detail"],
            )
            rejection_writer.writeheader()

            accepted_row_number = 0
            for input_row_number, row in enumerate(reader, 2):
                record_id = (row.get("id") or f"row_{input_row_number}").strip()
                try:
                    family = family_value(row)
                    sequence, length, pair_count, canonical_count = record_metrics(row)
                except ValueError as error:
                    if "non-AUCG" not in str(error) and "empty sequence" not in str(error):
                        raise ValueError(
                            f"invalid candidate row {input_row_number} ({record_id}): {error}"
                        ) from error
                    rejection_writer.writerow({
                        "id": record_id,
                        "family": row.get("family") or row.get("rfam_id") or "",
                        "reason": "non_aucg_sequence",
                        "detail": str(error),
                    })
                    continue

                failure = None
                values = reference.get(family)
                if values:
                    width = args.standard_deviations
                    length_min = max(0, values["length_mean"] - width * values["length_sd"])
                    length_max = values["length_mean"] + width * values["length_sd"]
                    pair_min = max(0, values["pair_mean"] - width * values["pair_sd"])
                    canonical_min = max(
                        0, values["canonical_mean"] - width * values["canonical_sd"]
                    )

                    if length < length_min or length > length_max:
                        failure = ("length_outlier", f"{length} not in [{length_min:.3f}, {length_max:.3f}]")
                    elif pair_count < pair_min:
                        failure = ("low_pair_count", f"{pair_count} < {pair_min:.3f}")
                    elif canonical_count < canonical_min:
                        failure = (
                            "low_canonical_pair_count",
                            f"{canonical_count} < {canonical_min:.3f}",
                        )

                if failure:
                    rejection_writer.writerow({
                        "id": record_id,
                        "family": family,
                        "reason": failure[0],
                        "detail": failure[1],
                    })
                    continue

                accepted_row_number += 1
                output = dict(row)
                output.update({
                    "qc_row_number": accepted_row_number,
                    "qc_sequence_length": length,
                    "qc_pair_count": pair_count,
                    "qc_canonical_pair_count": canonical_count,
                })
                preliminary_writer.writerow(output)

                interval = interval_fields(row)
                if interval and not args.skip_overlap_resolution:
                    accession, strand, start, end = interval
                    database.execute(
                        """
                        INSERT INTO intervals
                        (row_number, accession, strand, start, end, evalue, record_id)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            accepted_row_number,
                            accession,
                            strand,
                            start,
                            end,
                            parse_optional_float(row.get("evalue")),
                            record_id,
                        ),
                    )
                if accepted_row_number % 100000 == 0:
                    database.commit()
            database.commit()

        overlap_rejections = (
            {} if args.skip_overlap_resolution else resolve_overlaps(database)
        )
        database.close()

        with (
            (temporary_path / "accepted_pre_overlap.csv").open(
                newline="", encoding="utf-8"
            ) as preliminary_handle,
            args.output.open("w", newline="", encoding="utf-8") as output_handle,
            args.rejections.open("a", newline="", encoding="utf-8") as reject_handle,
        ):
            reader = csv.DictReader(preliminary_handle)
            output_fields = [field for field in reader.fieldnames if field != "qc_row_number"]
            writer = csv.DictWriter(output_handle, fieldnames=output_fields)
            writer.writeheader()
            rejection_writer = csv.DictWriter(
                reject_handle,
                fieldnames=["id", "family", "reason", "detail"],
            )
            for row in reader:
                row_number = int(row.pop("qc_row_number"))
                winner = overlap_rejections.get(row_number)
                if winner:
                    rejection_writer.writerow({
                        "id": row.get("id", ""),
                        "family": row.get("family") or row.get("rfam_id") or "",
                        "reason": "overlapping_hit",
                        "detail": f"retained {winner}",
                    })
                else:
                    writer.writerow(row)


if __name__ == "__main__":
    main()
