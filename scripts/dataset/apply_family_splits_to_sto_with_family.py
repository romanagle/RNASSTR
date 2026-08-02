#!/usr/bin/env python3

"""Convert partitioned Rfam Stockholm alignments to RNASSTR CSV files."""

import argparse
import csv
import json
import re
from pathlib import Path
from collections import defaultdict, Counter


GAP_CHARS = set(".-~_")
VALID_BASES = set("ACGUTNacgutn")


OPEN_TO_CLOSE = {
    "(": ")",
    "[": "]",
    "{": "}",
    "<": ">",
    "A": "a",
    "B": "b",
    "C": "c",
    "D": "d",
    "E": "e",
    "F": "f",
    "G": "g",
    "H": "h",
    "I": "i",
    "J": "j",
    "K": "k",
    "L": "l",
    "M": "m",
    "N": "n",
    "O": "o",
    "P": "p",
    "Q": "q",
    "R": "r",
    "S": "s",
    "T": "t",
    "U": "u",
    "V": "v",
    "W": "w",
    "X": "x",
    "Y": "y",
    "Z": "z",
}

CLOSE_TO_OPEN = {v: k for k, v in OPEN_TO_CLOSE.items()}


def normalize_family_name(text):
    """
    Extract RFxxxxx from a line, filename, or path.

    Examples:
        RF00005 -> RF00005
        RF00005.sto -> RF00005
        ./foo/RF02543_part01.sto -> RF02543
    """
    text = str(text).strip()
    if not text or text.startswith("#"):
        return None

    m = re.search(r"(RF\d{5})", text)
    if m:
        return m.group(1)

    return text


def read_family_list(path):
    families = set()
    path = Path(path)

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fam = normalize_family_name(line)
            if fam:
                families.add(fam)

    return families


def family_from_sto_path(path):
    return normalize_family_name(Path(path).name)


def clean_sequence(aligned_seq):
    """
    Remove alignment gaps and convert T to U.
    Keep ambiguous RNA/DNA letters as uppercase.
    """
    out = []
    for c in aligned_seq:
        if c in GAP_CHARS:
            continue

        c = c.upper()
        if c == "T":
            c = "U"

        out.append(c)

    return "".join(out)


def normalize_wuss_char(c):
    """
    Convert WUSS-ish structural chars to dot-bracket-like symbols.

    Keeps bracket characters that can be paired.
    Converts common unpaired/annotation chars to dot.
    """
    if c in OPEN_TO_CLOSE or c in CLOSE_TO_OPEN:
        return c

    if c in ".:,;_-~":
        return "."

    return "."


def parse_alignment_pairs(ss_cons):
    """
    Parse base pairs from the full alignment-level SS_cons string.

    Returns:
        dict mapping alignment_column_i -> alignment_column_j

    Handles multiple bracket types and simple pseudoknot letters.
    Unmatched brackets are ignored.
    """
    stacks = defaultdict(list)
    pairmap = {}

    for i, raw_c in enumerate(ss_cons):
        c = normalize_wuss_char(raw_c)

        if c in OPEN_TO_CLOSE:
            stacks[c].append(i)

        elif c in CLOSE_TO_OPEN:
            opener = CLOSE_TO_OPEN[c]
            if stacks[opener]:
                j = stacks[opener].pop()
                pairmap[i] = j
                pairmap[j] = i

    return pairmap


def project_structure_to_sequence(aligned_seq, ss_cons, pair_index_base=0):
    """
    Project alignment consensus structure onto an individual ungapped sequence.

    If one side of a consensus base pair is gapped in this sequence, that
    position is treated as unpaired in the projected structure.
    """
    if not ss_cons:
        sequence = clean_sequence(aligned_seq)
        structure = "." * len(sequence)
        return sequence, structure, []

    if len(ss_cons) < len(aligned_seq):
        ss_cons = ss_cons + "." * (len(aligned_seq) - len(ss_cons))
    elif len(ss_cons) > len(aligned_seq):
        ss_cons = ss_cons[: len(aligned_seq)]

    pairmap = parse_alignment_pairs(ss_cons)

    aln_to_seq = {}
    seq_chars = []

    seq_i = 0
    for aln_i, base in enumerate(aligned_seq):
        if base in GAP_CHARS:
            continue

        b = base.upper()
        if b == "T":
            b = "U"

        seq_chars.append(b)
        aln_to_seq[aln_i] = seq_i
        seq_i += 1

    sequence = "".join(seq_chars)
    struct_chars = ["." for _ in sequence]
    base_pairs = []

    for aln_i, aln_j in pairmap.items():
        if aln_i >= aln_j:
            continue

        if aln_i not in aln_to_seq or aln_j not in aln_to_seq:
            continue

        i = aln_to_seq[aln_i]
        j = aln_to_seq[aln_j]

        struct_chars[i] = "("
        struct_chars[j] = ")"

        if pair_index_base == 1:
            base_pairs.append([i + 1, j + 1])
        else:
            base_pairs.append([i, j])

    return sequence, "".join(struct_chars), base_pairs


def parse_stockholm(path):
    """
    Parse a possibly wrapped Stockholm file.

    Returns:
        records: dict seq_id -> aligned_sequence
        gc: dict annotation_name -> annotation_string

    Specifically captures:
        #=GC SS_cons
        #=GC RF
        and any other #=GC annotation encountered.

    Sequence lines are accumulated across blocks.
    """
    records = defaultdict(list)
    gc = defaultdict(list)

    with Path(path).open(errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if not line.strip():
                continue

            if line.startswith("# STOCKHOLM"):
                continue

            if line.strip() == "//":
                continue

            if line.startswith("#=GC"):
                parts = line.split(maxsplit=2)
                if len(parts) >= 3:
                    _, tag, chunk = parts
                    gc[tag].append(chunk.strip())
                continue

            if line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            seq_id = parts[0]
            chunk = parts[1]
            records[seq_id].append(chunk)

    records = {seq_id: "".join(chunks) for seq_id, chunks in records.items()}
    gc = {tag: "".join(chunks) for tag, chunks in gc.items()}

    return records, gc


def make_output_id(family, seq_id, id_mode):
    if id_mode == "seqname":
        return seq_id
    if id_mode == "family_seq":
        return f"{family}|{seq_id}"
    raise ValueError(f"Unknown id_mode: {id_mode}")


def discover_sto_files(sto_dir, recursive=True):
    sto_dir = Path(sto_dir)
    patterns = ["*.sto", "*.stk", "*.stockholm"]

    files = []
    for pat in patterns:
        if recursive:
            files.extend(sto_dir.rglob(pat))
        else:
            files.extend(sto_dir.glob(pat))

    return sorted(set(files))


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["id", "sequence", "structure", "base_pairs", "len", "family"]

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Apply family-level train/val/test splits to Stockholm files and write CSVs."
    )

    parser.add_argument("--sto-dir", required=True, help="Directory containing .sto files.")
    parser.add_argument("--train-families", required=True, help="Text file of train RFxxxxx families.")
    parser.add_argument("--val-families", required=True, help="Text file of val RFxxxxx families.")
    parser.add_argument("--test-families", required=True, help="Text file of test RFxxxxx families.")
    parser.add_argument("--outdir", required=True, help="Output directory for train.csv, val.csv, test.csv.")

    parser.add_argument(
        "--nonrecursive",
        action="store_true",
        help="Only scan --sto-dir directly; do not recurse into subdirectories.",
    )

    parser.add_argument(
        "--pair-index-base",
        type=int,
        choices=[0, 1],
        default=0,
        help="Base-pair coordinate indexing. Default: 0-based.",
    )

    parser.add_argument(
        "--id-mode",
        choices=["seqname", "family_seq"],
        default="family_seq",
        help="How to name output records. Default: family_seq.",
    )

    parser.add_argument(
        "--debug-family",
        default=None,
        help="Optional RFxxxxx family to print extra debugging info for.",
    )

    args = parser.parse_args()

    sto_dir = Path(args.sto_dir)
    outdir = Path(args.outdir)

    train_fams = read_family_list(args.train_families)
    val_fams = read_family_list(args.val_families)
    test_fams = read_family_list(args.test_families)

    split_for_family = {}

    overlaps = []
    for fam in train_fams:
        split_for_family[fam] = "train"
    for fam in val_fams:
        if fam in split_for_family:
            overlaps.append(fam)
        split_for_family[fam] = "val"
    for fam in test_fams:
        if fam in split_for_family:
            overlaps.append(fam)
        split_for_family[fam] = "test"

    if overlaps:
        print("WARNING: Some families appear in multiple split files.")
        print("Later split files override earlier ones in this order: train, val, test.")
        print("Overlapping families:", ", ".join(sorted(set(overlaps))[:50]))
        if len(set(overlaps)) > 50:
            print(f"... plus {len(set(overlaps)) - 50} more")

    sto_files = discover_sto_files(sto_dir, recursive=not args.nonrecursive)

    rows_by_split = {
        "train": [],
        "val": [],
        "test": [],
    }

    file_counts_by_split = Counter()
    seq_counts_by_family = Counter()
    skipped_no_family = []
    skipped_not_in_split = []
    skipped_parse_empty = []
    skipped_empty_records = 0

    duplicate_ids = Counter()

    for sto_path in sto_files:
        family = family_from_sto_path(sto_path)

        if family is None:
            skipped_no_family.append(str(sto_path))
            continue

        split = split_for_family.get(family)

        if split is None:
            skipped_not_in_split.append((family, str(sto_path)))
            continue

        try:
            records, gc = parse_stockholm(sto_path)
        except Exception as e:
            skipped_parse_empty.append((family, str(sto_path), f"PARSE_ERROR: {e}"))
            continue

        if not records:
            skipped_parse_empty.append((family, str(sto_path), "NO_SEQUENCE_RECORDS"))
            continue

        ss_cons = gc.get("SS_cons", "")

        if args.debug_family and family == args.debug_family:
            print("\nDEBUG FAMILY:", family)
            print("  file:", sto_path)
            print("  n_records:", len(records))
            print("  has_SS_cons:", bool(ss_cons))
            print("  SS_cons_len:", len(ss_cons))
            first_id = next(iter(records))
            print("  first_seq_id:", first_id)
            print("  first_aligned_len:", len(records[first_id]))
            print("  first_aligned_head:", records[first_id][:80])
            if ss_cons:
                print("  SS_cons_head:", ss_cons[:80])

        file_counts_by_split[split] += 1

        for seq_id, aligned_seq in records.items():
            sequence, structure, base_pairs = project_structure_to_sequence(
                aligned_seq,
                ss_cons,
                pair_index_base=args.pair_index_base,
            )

            if not sequence:
                skipped_empty_records += 1
                continue

            out_id = make_output_id(family, seq_id, args.id_mode)
            duplicate_ids[out_id] += 1

            if duplicate_ids[out_id] > 1:
                out_id = f"{out_id}|dup{duplicate_ids[out_id]}"

            row = {
                "id": out_id,
                "sequence": sequence,
                "structure": structure,
                "base_pairs": json.dumps(base_pairs, separators=(",", ":")),
                "len": len(sequence),
                "family": family,
            }

            rows_by_split[split].append(row)
            seq_counts_by_family[family] += 1

    outdir.mkdir(parents=True, exist_ok=True)

    write_csv(outdir / "train.csv", rows_by_split["train"])
    write_csv(outdir / "val.csv", rows_by_split["val"])
    write_csv(outdir / "test.csv", rows_by_split["test"])

    summary_path = outdir / "csv_split_summary.tsv"
    with summary_path.open("w") as handle:
        handle.write("split\tn_families_in_split_file\tn_sto_files_used\tn_sequences_written\n")
        handle.write(f"train\t{len(train_fams)}\t{file_counts_by_split['train']}\t{len(rows_by_split['train'])}\n")
        handle.write(f"val\t{len(val_fams)}\t{file_counts_by_split['val']}\t{len(rows_by_split['val'])}\n")
        handle.write(f"test\t{len(test_fams)}\t{file_counts_by_split['test']}\t{len(rows_by_split['test'])}\n")

    family_counts_path = outdir / "family_sequence_counts.tsv"
    with family_counts_path.open("w") as handle:
        handle.write("family\tsplit\tn_sequences\n")
        for fam in sorted(seq_counts_by_family):
            handle.write(f"{fam}\t{split_for_family.get(fam, 'NA')}\t{seq_counts_by_family[fam]}\n")

    skipped_path = outdir / "skipped_files.tsv"
    with skipped_path.open("w") as handle:
        handle.write("reason\tfamily\tpath\tdetail\n")

        for path in skipped_no_family:
            handle.write(f"NO_RFAMILY_IN_FILENAME\tNA\t{path}\tNA\n")

        for fam, path in skipped_not_in_split:
            handle.write(f"FAMILY_NOT_IN_SPLIT_FILES\t{fam}\t{path}\tNA\n")

        for fam, path, detail in skipped_parse_empty:
            handle.write(f"PARSE_EMPTY_OR_ERROR\t{fam}\t{path}\t{detail}\n")

    print("\nCSV files written:")
    print(f"  train: {len(rows_by_split['train'])} sequences -> {outdir / 'train.csv'}")
    print(f"  val:   {len(rows_by_split['val'])} sequences -> {outdir / 'val.csv'}")
    print(f"  test:  {len(rows_by_split['test'])} sequences -> {outdir / 'test.csv'}")

    print("\nDiagnostics:")
    print(f"  STO files discovered: {len(sto_files)}")
    print(f"  STO files used:")
    print(f"    train: {file_counts_by_split['train']}")
    print(f"    val:   {file_counts_by_split['val']}")
    print(f"    test:  {file_counts_by_split['test']}")
    print(f"  Files skipped because no RFxxxxx was found in filename: {len(skipped_no_family)}")
    print(f"  Files skipped because family was not in split files: {len(skipped_not_in_split)}")
    print(f"  Files skipped because no records were parsed / parse error: {len(skipped_parse_empty)}")
    print(f"  Empty/all-gap sequence records skipped: {skipped_empty_records}")

    print("\nExtra files written:")
    print(f"  {summary_path}")
    print(f"  {family_counts_path}")
    print(f"  {skipped_path}")

    if len(rows_by_split["train"]) == 0 and len(rows_by_split["val"]) == 0 and len(rows_by_split["test"]) == 0:
        print("\nWARNING: zero total sequences were written.")
        print("Most likely causes:")
        print("  1. The .sto filenames do not contain RFxxxxx family IDs.")
        print("  2. The family split files do not contain IDs matching the .sto filenames.")
        print("  3. The files are not under --sto-dir or have a nonstandard extension.")
        print("Check skipped_files.tsv first.")


if __name__ == "__main__":
    main()
