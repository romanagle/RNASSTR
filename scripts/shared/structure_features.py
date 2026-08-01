# shared/structure_features.py

"""
RNA secondary structure feature extraction utilities.

Used throughout the RNASSTR Phase 1 analysis pipeline.

Implements:
- GC content
- fraction paired
- stem count
- base pair extraction
- pseudoknot detection
- sequence length
- structural complexity metrics

Supports:
- standard dot-bracket
- pseudoknot brackets

References:
ViennaRNA dot-bracket notation:
https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/io/rna_structures.html
 [oai_citation:0‡tbi.univie.ac.at](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/io/rna_structures.html?utm_source=chatgpt.com)
"""

import re
import numpy as np


# ============================================================
# Valid bracket pairs
# ============================================================

OPEN_TO_CLOSE = {
    "(": ")",
    "[": "]",
    "{": "}",
    "<": ">"
}

CLOSE_TO_OPEN = {
    ")": "(",
    "]": "[",
    "}": "{",
    ">": "<"
}

ALL_OPEN = set(OPEN_TO_CLOSE.keys())
ALL_CLOSE = set(CLOSE_TO_OPEN.keys())


# ============================================================
# Sequence metrics
# ============================================================

def sequence_length(sequence):
    """
    Return sequence length.
    """

    return len(sequence)


def gc_content(sequence):
    """
    Fraction GC content.

    Returns value from 0-1.
    """

    sequence = sequence.upper()

    gc = sequence.count("G") + sequence.count("C")

    if len(sequence) == 0:
        return np.nan

    return gc / len(sequence)


# ============================================================
# Dot-bracket parsing
# ============================================================

def dotbracket_to_pairs(structure):
    """
    Convert dot-bracket notation to base pair list.

    Supports pseudoknots:
    (), [], {}, <>

    Returns
    -------
    list of tuples:
        [(i, j), ...]
    """

    stacks = {
        "(": [],
        "[": [],
        "{": [],
        "<": []
    }

    pairs = []

    for i, char in enumerate(structure):

        # opening bracket
        if char in ALL_OPEN:

            stacks[char].append(i)

        # closing bracket
        elif char in ALL_CLOSE:

            opener = CLOSE_TO_OPEN[char]

            if len(stacks[opener]) == 0:
                continue

            j = stacks[opener].pop()

            pairs.append((j, i))

    pairs = sorted(pairs)

    return pairs


# ============================================================
# Fraction paired
# ============================================================

def fraction_paired(structure):
    """
    Fraction paired nucleotides.

    Returns value from 0-1.
    """

    paired_chars = set("()[]{}<>")

    paired = sum(c in paired_chars for c in structure)

    if len(structure) == 0:
        return np.nan

    return paired / len(structure)


# ============================================================
# Stem counting
# ============================================================

def count_stems(structure):
    """
    Count stems/helices in dot-bracket structure.

    A new stem begins whenever:
    current pair is not contiguous with previous pair.

    Example:
        (((...))) -> 1 stem

        (((...)))...(((...)))
        -> 2 stems
    """

    pairs = dotbracket_to_pairs(structure)

    if len(pairs) == 0:
        return 0

    stems = 1

    prev_i, prev_j = pairs[0]

    for i, j in pairs[1:]:

        contiguous = (
            i == prev_i + 1
            and
            j == prev_j - 1
        )

        if not contiguous:
            stems += 1

        prev_i, prev_j = i, j

    return stems


# ============================================================
# Pseudoknot detection
# ============================================================

def contains_pseudoknot(structure):
    """
    Detect whether structure contains PK notation.
    """

    pk_chars = set("[]{}<>")

    return any(c in pk_chars for c in structure)


# ============================================================
# Structural density
# ============================================================

def basepair_density(structure):
    """
    Number of basepairs normalized by length.
    """

    pairs = dotbracket_to_pairs(structure)

    if len(structure) == 0:
        return np.nan

    return len(pairs) / len(structure)


# ============================================================
# Structural feature vector
# ============================================================

def compute_structure_features(sequence, structure):
    """
    Compute unified feature dictionary.

    Returns
    -------
    dict
    """

    return {

        "length":
            sequence_length(sequence),

        "gc_content":
            gc_content(sequence),

        "fraction_paired":
            fraction_paired(structure),

        "stem_count":
            count_stems(structure),

        "basepair_density":
            basepair_density(structure),

        "num_basepairs":
            len(dotbracket_to_pairs(structure)),

        "has_pseudoknot":
            contains_pseudoknot(structure),
    }


# ============================================================
# Main testing
# ============================================================

if __name__ == "__main__":

    seq = "GGGAAAUCCCGGGAAA"

    struct = "(((....)))(((....)))"

    print(seq)
    print(struct)

    print()

    print(dotbracket_to_pairs(struct))

    print()

    feats = compute_structure_features(seq, struct)

    for k, v in feats.items():
        print(f"{k}: {v}")
