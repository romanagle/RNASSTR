#!/usr/bin/env python3

"""Create structural- and clan-aware Rfam family partitions."""

import argparse
import gzip
import math
import random
import re
import shlex
import subprocess
from pathlib import Path
from collections import defaultdict


RF_RE = re.compile(r"RF\d{5}")
GAP_CHARS = set(".-~_")


class DSU:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x):
        self.add(x)
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        self.add(a)
        self.add(b)

        ra = self.find(a)
        rb = self.find(b)

        if ra == rb:
            return False

        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1

        return True


def log(msg):
    print(msg, flush=True)


def extract_rf_ids(text):
    return RF_RE.findall(str(text))


def family_from_path(path):
    m = RF_RE.search(Path(path).name)
    if not m:
        return None
    return m.group(0)


def open_text(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return path.open("r", errors="replace")


def discover_sto_files(sto_dir, recursive=True):
    sto_dir = Path(sto_dir)
    patterns = [
        "*.sto",
        "*.stk",
        "*.stockholm",
        "*.sto.gz",
        "*.stk.gz",
        "*.stockholm.gz",
    ]

    files = []
    for pat in patterns:
        if recursive:
            files.extend(sto_dir.rglob(pat))
        else:
            files.extend(sto_dir.glob(pat))

    return sorted(set(files))


def clean_seq(aligned):
    out = []

    for c in aligned:
        if c in GAP_CHARS:
            continue

        c = c.upper()
        if c == "T":
            c = "U"

        out.append(c)

    return "".join(out)


def parse_stockholm_sequences(path):
    """
    Parse wrapped Stockholm sequence records.

    Returns:
        dict seq_id -> ungapped RNA sequence
    """
    aligned_chunks = defaultdict(list)

    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")

            if not line.strip():
                continue

            if line.strip() == "//":
                continue

            if line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            seq_id = parts[0]
            seq_chunk = parts[1]
            aligned_chunks[seq_id].append(seq_chunk)

    seqs = {}

    for seq_id, chunks in aligned_chunks.items():
        seq = clean_seq("".join(chunks))
        if seq:
            seqs[seq_id] = seq

    return seqs


def load_family_sequences(sto_dir, recursive=True):
    """
    Load all family sequence IDs and ungapped sequences from STO files.

    Families are inferred from RFxxxxx in filenames.
    Multiple files belonging to same family, e.g. RF02543_part01.sto, are merged.
    """
    sto_files = discover_sto_files(sto_dir, recursive=recursive)

    family_to_records = defaultdict(list)
    skipped_no_family = []
    skipped_empty = []

    for path in sto_files:
        fam = family_from_path(path)

        if fam is None:
            skipped_no_family.append(str(path))
            continue

        seqs = parse_stockholm_sequences(path)

        if not seqs:
            skipped_empty.append(str(path))
            continue

        for seq_id, seq in seqs.items():
            family_to_records[fam].append((seq_id, seq, str(path)))

    return family_to_records, sto_files, skipped_no_family, skipped_empty


def write_sampled_fasta(family_to_records, out_fasta, max_seqs_per_family, seed):
    """
    Write sampled sequences to FASTA.

    FASTA headers:
        >RFxxxxx|sample_index|original_seq_id
    """
    rng = random.Random(seed)

    total_written = 0
    sampled_counts = {}

    with Path(out_fasta).open("w") as out:
        for fam in sorted(family_to_records):
            records = list(family_to_records[fam])

            if max_seqs_per_family is not None and len(records) > max_seqs_per_family:
                records = rng.sample(records, max_seqs_per_family)

            sampled_counts[fam] = len(records)

            for i, (seq_id, seq, _src_path) in enumerate(records, start=1):
                safe_seq_id = re.sub(r"\s+", "_", seq_id)
                header = f"{fam}|sample{i}|{safe_seq_id}"

                out.write(f">{header}\n")

                for j in range(0, len(seq), 80):
                    out.write(seq[j:j + 80] + "\n")

                total_written += 1

    return sampled_counts, total_written


def write_family_counts(family_to_records, sampled_counts, out_path):
    with Path(out_path).open("w") as out:
        out.write("family\tn_sequences_total\tn_sequences_sampled\n")

        for fam in sorted(family_to_records):
            out.write(
                f"{fam}\t{len(family_to_records[fam])}\t"
                f"{sampled_counts.get(fam, 0)}\n"
            )


def ensure_cm_indexes(rfam_cm):
    rfam_cm = Path(rfam_cm)
    suffixes = [".i1f", ".i1i", ".i1m", ".i1p"]

    missing = []

    for suffix in suffixes:
        if not Path(str(rfam_cm) + suffix).exists():
            missing.append(str(rfam_cm) + suffix)

    return missing


def run_cmscan(
    cmscan_bin,
    rfam_cm,
    fasta,
    tblout,
    logfile,
    cpu,
    evalue,
    nohmmonly=True,
    extra_args=None,
):
    """
    Run Infernal cmscan.

    Note:
      cmscan supports --tblout.
      It does NOT support --domtblout.
    """
    cmd = [
        cmscan_bin,
        "--cpu",
        str(cpu),
        "-E",
        str(evalue),
        "--tblout",
        str(tblout),
    ]

    if nohmmonly:
        cmd.append("--nohmmonly")

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend([str(rfam_cm), str(fasta)])

    log("\nRunning cmscan:")
    log("  " + " ".join(shlex.quote(x) for x in cmd))

    with Path(logfile).open("w") as log_handle:
        proc = subprocess.run(
            cmd,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )

    if proc.returncode != 0:
        raise RuntimeError(
            f"cmscan failed with return code {proc.returncode}. "
            f"Check log file: {logfile}"
        )


def parse_float_safe(x):
    try:
        return float(x)
    except Exception:
        return None


def source_family_from_query_name(qname):
    fams = extract_rf_ids(qname)
    if fams:
        return fams[0]
    return None


def target_family_from_tbl_tokens(tokens):
    """
    Infer target RF family from cmscan --tblout tokens.

    Usually:
      target name, accession, query name, accession, ...
    Either target name or accession may include RFxxxxx.
    """
    if len(tokens) < 2:
        return None

    for idx in [0, 1]:
        fams = extract_rf_ids(tokens[idx])
        if fams:
            return fams[0]

    fams = extract_rf_ids(" ".join(tokens))
    if fams:
        return fams[0]

    return None


def query_name_from_tbl_tokens(tokens):
    """
    Infer query name from cmscan --tblout tokens.

    For cmscan tblout, query name is usually token 2.
    """
    if len(tokens) >= 3:
        return tokens[2]
    return None


def evalue_from_tbl_tokens(tokens):
    """
    Infer full-sequence E-value from cmscan --tblout.

    Typical cmscan --tblout has E-value around token 16,
    but this tries a few likely positions and then falls back.
    """
    likely_indices = [16, 15, 17, 18]

    for idx in likely_indices:
        if idx < len(tokens):
            val = parse_float_safe(tokens[idx])
            if val is not None:
                return val

    for tok in tokens[10:]:
        val = parse_float_safe(tok)
        if val is not None and val >= 0:
            return val

    return None


def parse_cmscan_tblout(tblout, evalue_threshold, valid_families):
    """
    Parse cmscan --tblout and return structural cross-family edges.

    Returns:
      edges: dict (source_family, target_family) -> best_evalue
      hit_rows: list of rows for writing diagnostics
    """
    edges = {}
    hit_rows = []

    with Path(tblout).open(errors="replace") as handle:
        for raw in handle:
            line = raw.strip()

            if not line or line.startswith("#"):
                continue

            tokens = line.split()

            if len(tokens) < 4:
                continue

            target = target_family_from_tbl_tokens(tokens)
            qname = query_name_from_tbl_tokens(tokens)
            source = source_family_from_query_name(qname) if qname else None
            evalue = evalue_from_tbl_tokens(tokens)

            if source is None or target is None:
                continue

            if source not in valid_families or target not in valid_families:
                continue

            if evalue is None:
                continue

            if evalue > evalue_threshold:
                continue

            hit_rows.append(
                {
                    "source_family": source,
                    "target_family": target,
                    "evalue": evalue,
                    "query_name": qname,
                    "target_token": tokens[0],
                    "raw": line,
                }
            )

            if source == target:
                continue

            a, b = sorted([source, target])
            key = (a, b)

            if key not in edges or evalue < edges[key]:
                edges[key] = evalue

    return edges, hit_rows


def parse_clanin(clanin, valid_families):
    """
    Parse Rfam.clanin loosely.

    Any non-comment line with two or more RFxxxxx IDs creates a clan group.
    """
    clan_groups = []

    if not clanin:
        return clan_groups

    with Path(clanin).open(errors="replace") as handle:
        for line_i, raw in enumerate(handle, start=1):
            line = raw.strip()

            if not line or line.startswith("#"):
                continue

            fams = sorted(set(f for f in extract_rf_ids(line) if f in valid_families))

            if len(fams) >= 2:
                clan_groups.append((f"clan_line_{line_i}", fams, raw.rstrip("\n")))

    return clan_groups


def build_components(valid_families, structural_edges, clan_groups):
    dsu = DSU()

    for fam in valid_families:
        dsu.add(fam)

    for a, b in structural_edges:
        dsu.union(a, b)

    for _clan_id, fams, _raw in clan_groups:
        first = fams[0]

        for fam in fams[1:]:
            dsu.union(first, fam)

    root_to_fams = defaultdict(list)

    for fam in valid_families:
        root_to_fams[dsu.find(fam)].append(fam)

    components = []

    for fams in root_to_fams.values():
        components.append(sorted(fams))

    components.sort(key=lambda xs: (len(xs), xs[0]))

    return components


def component_records(components, family_counts):
    records = []

    for i, fams in enumerate(components, start=1):
        nseq = sum(family_counts[f] for f in fams)

        records.append(
            {
                "component_id": f"C{i:05d}",
                "families": fams,
                "n_families": len(fams),
                "n_sequences": nseq,
            }
        )

    records.sort(key=lambda c: (c["n_sequences"], c["n_families"], c["families"][0]))

    for i, rec in enumerate(records, start=1):
        rec["component_id"] = f"C{i:05d}"

    return records


def split_size(comps):
    return sum(c["n_sequences"] for c in comps)


def split_nfams(comps):
    return sum(c["n_families"] for c in comps)


def choose_split_for_diversity(val_comps, test_comps):
    """
    Choose whether next small component goes to val or test.

    Balance by:
      1. family count
      2. component count
      3. sequence count
    """
    val_nfam = split_nfams(val_comps)
    test_nfam = split_nfams(test_comps)

    if val_nfam < test_nfam:
        return "val"
    if test_nfam < val_nfam:
        return "test"

    if len(val_comps) < len(test_comps):
        return "val"
    if len(test_comps) < len(val_comps):
        return "test"

    if split_size(val_comps) <= split_size(test_comps):
        return "val"

    return "test"


def is_large_component(comp, total_sequences, max_holdout_component_frac):
    if total_sequences <= 0:
        return False

    return (comp["n_sequences"] / total_sequences) > max_holdout_component_frac


def assign_components(
    components,
    train_frac,
    val_frac,
    test_frac,
    min_val_families,
    min_test_families,
    min_val_components,
    min_test_components,
    max_holdout_component_frac,
):
    """
    Diversity-first split assignment.

    Objective order:
      1. Keep structural/clan components intact.
      2. Put many small independent components/families into val/test.
      3. Avoid placing giant components in val/test.
      4. Only then try to get near target sequence fractions.
    """
    total = sum(c["n_sequences"] for c in components)

    val_target = total * val_frac
    test_target = total * test_frac

    small_components = []
    large_components = []

    for comp in components:
        if is_large_component(comp, total, max_holdout_component_frac):
            large_components.append(comp)
        else:
            small_components.append(comp)

    small_components = sorted(
        small_components,
        key=lambda c: (c["n_sequences"], c["n_families"], c["families"][0]),
    )

    large_components = sorted(
        large_components,
        key=lambda c: (c["n_sequences"], c["n_families"], c["families"][0]),
        reverse=True,
    )

    val_comps = []
    test_comps = []
    used_ids = set()

    for comp in small_components:
        val_needs_diversity = (
            split_nfams(val_comps) < min_val_families
            or len(val_comps) < min_val_components
        )

        test_needs_diversity = (
            split_nfams(test_comps) < min_test_families
            or len(test_comps) < min_test_components
        )

        if not val_needs_diversity and not test_needs_diversity:
            break

        if val_needs_diversity and test_needs_diversity:
            dest = choose_split_for_diversity(val_comps, test_comps)
        elif val_needs_diversity:
            dest = "val"
        else:
            dest = "test"

        if dest == "val":
            val_comps.append(comp)
        else:
            test_comps.append(comp)

        used_ids.add(comp["component_id"])

    remaining_small = [
        c for c in small_components
        if c["component_id"] not in used_ids
    ]

    for comp in remaining_small:
        val_size = split_size(val_comps)
        test_size = split_size(test_comps)

        val_under = val_size < val_target
        test_under = test_size < test_target

        if not val_under and not test_under:
            break

        if val_under and test_under:
            val_deficit = val_target - val_size
            test_deficit = test_target - test_size

            if val_deficit >= test_deficit:
                val_comps.append(comp)
            else:
                test_comps.append(comp)

        elif val_under:
            val_comps.append(comp)

        elif test_under:
            test_comps.append(comp)

        used_ids.add(comp["component_id"])

    train_comps = [
        c for c in components
        if c["component_id"] not in used_ids
    ]

    return train_comps, val_comps, test_comps


def polish_by_single_component_moves(
    train_comps,
    val_comps,
    test_comps,
    total,
    val_frac,
    test_frac,
    max_holdout_component_frac,
):
    """
    Diversity-safe polishing.

    This never removes components from val/test.
    It only adds small train components to val/test if doing so improves
    sequence-balance.
    """
    val_target = total * val_frac
    test_target = total * test_frac

    train_comps.sort(
        key=lambda c: (c["n_sequences"], c["n_families"], c["families"][0])
    )

    changed = True

    while changed:
        changed = False

        val_size = split_size(val_comps)
        test_size = split_size(test_comps)

        val_err = abs(val_size - val_target)
        test_err = abs(test_size - test_target)

        best_move = None
        best_dest = None
        best_improvement = 0

        for comp in train_comps:
            if total > 0 and (comp["n_sequences"] / total) > max_holdout_component_frac:
                continue

            new_val_err = abs((val_size + comp["n_sequences"]) - val_target)
            val_improvement = val_err - new_val_err

            if val_improvement > best_improvement:
                best_improvement = val_improvement
                best_move = comp
                best_dest = "val"

            new_test_err = abs((test_size + comp["n_sequences"]) - test_target)
            test_improvement = test_err - new_test_err

            if test_improvement > best_improvement:
                best_improvement = test_improvement
                best_move = comp
                best_dest = "test"

        if best_move is not None and best_improvement > 0:
            train_comps.remove(best_move)

            if best_dest == "val":
                val_comps.append(best_move)
            else:
                test_comps.append(best_move)

            changed = True

    return train_comps, val_comps, test_comps


def comps_to_families(comps):
    fams = []

    for comp in comps:
        fams.extend(comp["families"])

    return sorted(fams)


def write_family_file(path, fams):
    with Path(path).open("w") as out:
        for fam in sorted(fams):
            out.write(f"{fam}\n")


def write_structural_hits(path, hit_rows):
    with Path(path).open("w") as out:
        out.write("source_family\ttarget_family\tevalue\tquery_name\ttarget_token\traw\n")

        for row in hit_rows:
            out.write(
                f"{row['source_family']}\t{row['target_family']}\t"
                f"{row['evalue']:.6g}\t{row['query_name']}\t"
                f"{row['target_token']}\t{row['raw']}\n"
            )


def write_structural_edges(path, edges):
    with Path(path).open("w") as out:
        out.write("family_a\tfamily_b\tbest_evalue\n")

        for (a, b), evalue in sorted(edges.items()):
            out.write(f"{a}\t{b}\t{evalue:.6g}\n")


def write_clan_edges(path, clan_groups):
    with Path(path).open("w") as out:
        out.write("clan_id\tfamilies\toriginal_line\n")

        for clan_id, fams, raw in clan_groups:
            out.write(f"{clan_id}\t{','.join(fams)}\t{raw}\n")


def write_component_summary(path, split_to_comps):
    with Path(path).open("w") as out:
        out.write("component_id\tsplit\tn_families\tn_sequences\tfamilies\n")

        for split in ["train", "val", "test"]:
            comps = split_to_comps[split]

            for comp in sorted(comps, key=lambda c: c["component_id"]):
                out.write(
                    f"{comp['component_id']}\t{split}\t"
                    f"{comp['n_families']}\t{comp['n_sequences']}\t"
                    f"{','.join(comp['families'])}\n"
                )


def write_split_summary(path, split_to_comps):
    total = sum(split_size(comps) for comps in split_to_comps.values())

    with Path(path).open("w") as out:
        out.write("split\tn_components\tn_families\tn_sequences\tsequence_fraction\n")

        for split in ["train", "val", "test"]:
            comps = split_to_comps[split]
            nseq = split_size(comps)
            nfam = split_nfams(comps)
            frac = nseq / total if total else 0

            out.write(
                f"{split}\t{len(comps)}\t{nfam}\t{nseq}\t{frac:.8f}\n"
            )


def write_family_assignments(path, split_to_comps, family_counts):
    fam_to_split = {}

    for split, comps in split_to_comps.items():
        for comp in comps:
            for fam in comp["families"]:
                fam_to_split[fam] = split

    with Path(path).open("w") as out:
        out.write("family\tsplit\tn_sequences\n")

        for fam in sorted(family_counts):
            out.write(f"{fam}\t{fam_to_split.get(fam, 'NA')}\t{family_counts[fam]}\n")


def write_clan_split_check(path, clan_groups, split_to_comps):
    fam_to_split = {}

    for split, comps in split_to_comps.items():
        for comp in comps:
            for fam in comp["families"]:
                fam_to_split[fam] = split

    with Path(path).open("w") as out:
        out.write("clan_id\tfamilies_present\tsplits_present\tstatus\toriginal_line\n")

        for clan_id, fams, raw in clan_groups:
            present = [f for f in fams if f in fam_to_split]
            splits = sorted(set(fam_to_split[f] for f in present))

            status = "OK" if len(splits) <= 1 else "VIOLATION"

            out.write(
                f"{clan_id}\t{','.join(present)}\t"
                f"{','.join(splits)}\t{status}\t{raw}\n"
            )


def write_structural_split_check(path, structural_edges, split_to_comps):
    fam_to_split = {}

    for split, comps in split_to_comps.items():
        for comp in comps:
            for fam in comp["families"]:
                fam_to_split[fam] = split

    with Path(path).open("w") as out:
        out.write("family_a\tfamily_b\tsplit_a\tsplit_b\tstatus\tbest_evalue\n")

        for (a, b), evalue in sorted(structural_edges.items()):
            split_a = fam_to_split.get(a, "NA")
            split_b = fam_to_split.get(b, "NA")
            status = "OK" if split_a == split_b else "VIOLATION"

            out.write(
                f"{a}\t{b}\t{split_a}\t{split_b}\t{status}\t{evalue:.6g}\n"
            )


def maybe_load_existing_edges(path, valid_families):
    """
    Load previously written structural_edges.tsv.

    Expected columns:
        family_a family_b best_evalue

    Also works loosely by extracting first two RF IDs per line.
    """
    edges = {}

    with Path(path).open(errors="replace") as handle:
        for raw in handle:
            line = raw.strip()

            if not line or line.startswith("#"):
                continue

            if line.lower().startswith("family_a"):
                continue

            fams = [f for f in extract_rf_ids(line) if f in valid_families]

            if len(fams) < 2:
                continue

            a, b = sorted(fams[:2])

            if a == b:
                continue

            parts = line.split()
            evalue = None

            for p in parts:
                val = parse_float_safe(p)
                if val is not None:
                    evalue = val
                    break

            if evalue is None:
                evalue = math.nan

            key = (a, b)

            if key not in edges:
                edges[key] = evalue
            elif not math.isnan(evalue) and (
                math.isnan(edges[key]) or evalue < edges[key]
            ):
                edges[key] = evalue

    return edges


def main():
    parser = argparse.ArgumentParser(
        description="Make diversity-first structural + clan-aware train/val/test splits."
    )

    parser.add_argument(
        "--sto-dir",
        required=True,
        help="Directory containing .sto files.",
    )

    parser.add_argument(
        "--rfam-cm",
        required=True,
        help="Rfam.cm covariance model file.",
    )

    parser.add_argument(
        "--clanin",
        required=True,
        help="Rfam.clanin file.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )

    parser.add_argument(
        "--max-seqs-per-family",
        type=int,
        default=10,
        help="Max sequences sampled per family for cmscan. Use 0 for all sequences. Default: 10.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed for per-family subsampling. Default: 1.",
    )

    parser.add_argument(
        "--evalue",
        type=float,
        default=0.01,
        help="E-value threshold for structural cross-family edges. Default: 0.01.",
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=64,
        help="Threads for cmscan. Default: 64.",
    )

    parser.add_argument(
        "--cmscan-bin",
        default="cmscan",
        help="cmscan executable. Default: cmscan.",
    )

    parser.add_argument(
        "--train-frac",
        type=float,
        default=0.90,
        help="Soft target train sequence fraction. Default: 0.90.",
    )

    parser.add_argument(
        "--val-frac",
        type=float,
        default=0.05,
        help="Soft target validation sequence fraction. Default: 0.05.",
    )

    parser.add_argument(
        "--test-frac",
        type=float,
        default=0.05,
        help="Soft target test sequence fraction. Default: 0.05.",
    )

    parser.add_argument(
        "--min-val-families",
        type=int,
        default=300,
        help="Minimum validation families if possible. Default: 300.",
    )

    parser.add_argument(
        "--min-test-families",
        type=int,
        default=300,
        help="Minimum test families if possible. Default: 300.",
    )

    parser.add_argument(
        "--min-val-components",
        type=int,
        default=150,
        help="Minimum validation structural/clan components if possible. Default: 150.",
    )

    parser.add_argument(
        "--min-test-components",
        type=int,
        default=150,
        help="Minimum test structural/clan components if possible. Default: 150.",
    )

    parser.add_argument(
        "--max-holdout-component-frac",
        type=float,
        default=0.01,
        help=(
            "Components larger than this total-sequence fraction are kept out "
            "of val/test during diversity filling. Default: 0.01."
        ),
    )

    parser.add_argument(
        "--nonrecursive",
        action="store_true",
        help="Only scan --sto-dir directly, not recursively.",
    )

    parser.add_argument(
        "--skip-cmscan",
        action="store_true",
        help="Skip cmscan and use --existing-edges instead.",
    )

    parser.add_argument(
        "--existing-edges",
        default=None,
        help="Existing structural_edges.tsv to reuse.",
    )

    parser.add_argument(
        "--hmm-only",
        action="store_true",
        help="Do not use cmscan --nohmmonly. Faster but less structure-sensitive.",
    )

    parser.add_argument(
        "--no-polish",
        action="store_true",
        help="Skip conservative post-hoc sequence-balance polishing.",
    )

    parser.add_argument(
        "--cmscan-extra-args",
        default="",
        help="Extra arguments passed to cmscan as a quoted string.",
    )

    args = parser.parse_args()

    frac_sum = args.train_frac + args.val_frac + args.test_frac

    if abs(frac_sum - 1.0) > 1e-8:
        raise ValueError(f"Fractions must sum to 1.0; got {frac_sum}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sampled_fasta = outdir / "sampled_sequences_for_cmscan.fa"
    family_counts_tsv = outdir / "family_counts.tsv"
    cmscan_tblout = outdir / "cmscan.tblout"
    cmscan_log = outdir / "cmscan.log"
    structural_hits_tsv = outdir / "structural_hits.tsv"
    structural_edges_tsv = outdir / "structural_edges.tsv"
    clan_edges_tsv = outdir / "clan_edges.tsv"
    component_summary_tsv = outdir / "component_summary.tsv"
    split_summary_tsv = outdir / "split_summary.tsv"
    family_assignments_tsv = outdir / "family_assignments.tsv"
    clan_check_tsv = outdir / "clan_split_check.tsv"
    structural_check_tsv = outdir / "structural_split_check.tsv"

    log("\nLoading Stockholm files...")

    family_to_records, sto_files, skipped_no_family, skipped_empty = load_family_sequences(
        args.sto_dir,
        recursive=not args.nonrecursive,
    )

    if not family_to_records:
        raise RuntimeError("No sequence records loaded from STO files.")

    family_counts = {
        fam: len(records)
        for fam, records in family_to_records.items()
    }

    valid_families = set(family_counts)

    log(f"  STO files discovered: {len(sto_files)}")
    log(f"  RF families loaded:    {len(family_to_records)}")
    log(f"  Total sequences:       {sum(family_counts.values())}")
    log(f"  Skipped no-family:     {len(skipped_no_family)}")
    log(f"  Skipped empty:         {len(skipped_empty)}")

    max_seqs = args.max_seqs_per_family

    if max_seqs is not None and max_seqs <= 0:
        max_seqs = None

    log("\nWriting sampled FASTA for cmscan...")

    sampled_counts, total_sampled = write_sampled_fasta(
        family_to_records=family_to_records,
        out_fasta=sampled_fasta,
        max_seqs_per_family=max_seqs,
        seed=args.seed,
    )

    write_family_counts(family_to_records, sampled_counts, family_counts_tsv)

    log(f"  Sampled FASTA: {sampled_fasta}")
    log(f"  Total sampled sequences: {total_sampled}")
    log(f"  Family counts: {family_counts_tsv}")

    structural_edges = {}
    hit_rows = []

    if args.skip_cmscan:
        if not args.existing_edges:
            raise ValueError("--skip-cmscan requires --existing-edges")

        log("\nLoading existing structural edges...")
        structural_edges = maybe_load_existing_edges(args.existing_edges, valid_families)
        log(f"  Loaded structural edges: {len(structural_edges)}")

    else:
        missing_indexes = ensure_cm_indexes(args.rfam_cm)

        if missing_indexes:
            log("\nWARNING: Missing cmpress index files for Rfam.cm:")

            for p in missing_indexes:
                log(f"  {p}")

            log("\nRun this before rerunning the script:")
            log(f"  cmpress {args.rfam_cm}")

            raise RuntimeError("Missing Rfam.cm index files.")

        extra_args = shlex.split(args.cmscan_extra_args) if args.cmscan_extra_args else []

        use_nohmmonly = not args.hmm_only

        run_cmscan(
            cmscan_bin=args.cmscan_bin,
            rfam_cm=args.rfam_cm,
            fasta=sampled_fasta,
            tblout=cmscan_tblout,
            logfile=cmscan_log,
            cpu=args.cpu,
            evalue=args.evalue,
            nohmmonly=use_nohmmonly,
            extra_args=extra_args,
        )

        log("\nParsing cmscan hits...")

        structural_edges, hit_rows = parse_cmscan_tblout(
            tblout=cmscan_tblout,
            evalue_threshold=args.evalue,
            valid_families=valid_families,
        )

        write_structural_hits(structural_hits_tsv, hit_rows)

    write_structural_edges(structural_edges_tsv, structural_edges)

    log(f"  Significant cross-family structural edges: {len(structural_edges)}")
    log(f"  Structural edges written: {structural_edges_tsv}")

    if hit_rows:
        log(f"  Structural hits written:  {structural_hits_tsv}")

    log("\nParsing Rfam clan constraints...")

    clan_groups = parse_clanin(args.clanin, valid_families)
    write_clan_edges(clan_edges_tsv, clan_groups)

    n_clan_family_memberships = sum(len(fams) for _cid, fams, _raw in clan_groups)

    log(f"  Clan groups with >=2 present families: {len(clan_groups)}")
    log(f"  Clan family memberships represented:  {n_clan_family_memberships}")
    log(f"  Clan edges/groups written: {clan_edges_tsv}")

    log("\nBuilding structural + clan connected components...")

    components_raw = build_components(
        valid_families=valid_families,
        structural_edges=structural_edges.keys(),
        clan_groups=clan_groups,
    )

    components = component_records(components_raw, family_counts)

    log(f"  Components: {len(components)}")
    log("  Largest components by sequence count:")

    for comp in sorted(components, key=lambda c: c["n_sequences"], reverse=True)[:10]:
        preview = ",".join(comp["families"][:8])

        if len(comp["families"]) > 8:
            preview += ",..."

        log(
            f"    {comp['component_id']}: "
            f"{comp['n_sequences']} seqs, "
            f"{comp['n_families']} families: {preview}"
        )

    log("\nAssigning components to train/val/test with diversity-first objective...")

    train_comps, val_comps, test_comps = assign_components(
        components=components,
        train_frac=args.train_frac,
        val_frac=args.val_frac,
        test_frac=args.test_frac,
        min_val_families=args.min_val_families,
        min_test_families=args.min_test_families,
        min_val_components=args.min_val_components,
        min_test_components=args.min_test_components,
        max_holdout_component_frac=args.max_holdout_component_frac,
    )

    total_sequences = sum(family_counts.values())

    if not args.no_polish:
        train_comps, val_comps, test_comps = polish_by_single_component_moves(
            train_comps=train_comps,
            val_comps=val_comps,
            test_comps=test_comps,
            total=total_sequences,
            val_frac=args.val_frac,
            test_frac=args.test_frac,
            max_holdout_component_frac=args.max_holdout_component_frac,
        )

    split_to_comps = {
        "train": train_comps,
        "val": val_comps,
        "test": test_comps,
    }

    train_fams = comps_to_families(train_comps)
    val_fams = comps_to_families(val_comps)
    test_fams = comps_to_families(test_comps)

    write_family_file(outdir / "train_families.txt", train_fams)
    write_family_file(outdir / "val_families.txt", val_fams)
    write_family_file(outdir / "test_families.txt", test_fams)

    write_component_summary(component_summary_tsv, split_to_comps)
    write_split_summary(split_summary_tsv, split_to_comps)
    write_family_assignments(family_assignments_tsv, split_to_comps, family_counts)
    write_clan_split_check(clan_check_tsv, clan_groups, split_to_comps)
    write_structural_split_check(structural_check_tsv, structural_edges, split_to_comps)

    log("\nSplit files written:")
    log(f"  {outdir / 'train_families.txt'}")
    log(f"  {outdir / 'val_families.txt'}")
    log(f"  {outdir / 'test_families.txt'}")

    log("\nFinal split summary:")

    for split in ["train", "val", "test"]:
        comps = split_to_comps[split]
        nseq = split_size(comps)
        nfam = split_nfams(comps)
        frac = nseq / total_sequences if total_sequences else 0

        log(
            f"  {split:5s}: "
            f"{len(comps):6d} components, "
            f"{nfam:6d} families, "
            f"{nseq:10d} sequences, "
            f"{frac:.4%}"
        )

    log("\nDiagnostics written:")
    log(f"  {family_counts_tsv}")
    log(f"  {sampled_fasta}")
    log(f"  {cmscan_tblout}")
    log(f"  {cmscan_log}")
    log(f"  {structural_edges_tsv}")
    log(f"  {clan_edges_tsv}")
    log(f"  {component_summary_tsv}")
    log(f"  {split_summary_tsv}")
    log(f"  {family_assignments_tsv}")
    log(f"  {clan_check_tsv}")
    log(f"  {structural_check_tsv}")

    log("\nUseful checks:")
    log(f"  cat {split_summary_tsv}")
    log(f"  grep VIOLATION {clan_check_tsv} | head")
    log(f"  grep VIOLATION {structural_check_tsv} | head")
    log(f"  grep -E 'RF00001|RF02547' {family_assignments_tsv}")
    log(f"  grep -E 'RF00001|RF02547' {component_summary_tsv}")


if __name__ == "__main__":
    main()
