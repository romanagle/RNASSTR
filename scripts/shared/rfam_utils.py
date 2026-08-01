# shared/rfam_utils.py

"""
Utilities for parsing:
- Rfam.cm
- Rfam.clanin

Extracts:
- RFAM accession
- short NAME
- CLEN
- NSEQ
- clan assignments

Used throughout the RNASSTR revision pipeline.
"""

from pathlib import Path
import pandas as pd


# ============================================================
# Rfam.cm parsing
# ============================================================

def parse_rfam_cm(cm_file):
    """
    Parse metadata from Rfam.cm.

    Extracts:
    - ACC
    - NAME
    - CLEN
    - NSEQ

    Returns
    -------
    pandas.DataFrame
    """

    cm_file = Path(cm_file)

    records = []

    current = {}

    with open(cm_file) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            # ------------------------------------------------
            # End of model block
            # ------------------------------------------------
            if line == "//":

                # only append COMPLETE records
                if "family" in current:

                    records.append(current.copy())

                current = {}

                continue

            # ------------------------------------------------
            # Parse fields
            # ------------------------------------------------

            if line.startswith("ACC"):

                fields = line.split(maxsplit=1)

                if len(fields) > 1:
                    current["family"] = fields[1].strip()

            elif line.startswith("NAME"):

                fields = line.split(maxsplit=1)

                if len(fields) > 1:
                    current["name"] = fields[1].strip()

            elif line.startswith("CLEN"):

                fields = line.split(maxsplit=1)

                if len(fields) > 1:

                    try:
                        current["clen"] = int(fields[1])

                    except:
                        pass

            elif line.startswith("NSEQ"):

                fields = line.split(maxsplit=1)

                if len(fields) > 1:

                    try:
                        current["nseq"] = int(fields[1])

                    except:
                        pass

    df = pd.DataFrame(records)

    # --------------------------------------------------------
    # Deduplicate families
    # --------------------------------------------------------

    # Prefer rows with non-null CLEN
    df["clen_missing"] = df["clen"].isna()

    df = (
        df
        .sort_values("clen_missing")
        .drop_duplicates(subset=["family"], keep="first")
        .drop(columns=["clen_missing"])
    )

    df = df.reset_index(drop=True)

    print(f"\nParsed {len(df)} unique RFAM families from {cm_file}")

    print(df.head())

    return df


# ============================================================
# Clan parsing
# ============================================================

def parse_rfam_clans(clan_file):
    """
    Parse Rfam.clanin.

    Returns
    -------
    pandas.DataFrame

    Columns:
    - family
    - clan
    """

    clan_file = Path(clan_file)

    records = []

    with open(clan_file) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 2:
                continue

            clan = fields[0]
            family = fields[1]

            records.append({
                "family": family,
                "clan": clan
            })

    df = pd.DataFrame(records)

    print(f"\nParsed {len(df)} clan mappings from {clan_file}")

    return df


# ============================================================
# Combined metadata table
# ============================================================

def build_rfam_metadata(cm_file, clan_file=None):
    """
    Build combined RFAM metadata table.
    """

    rfam_df = parse_rfam_cm(cm_file)

    if rfam_df.empty:

        raise ValueError(
            "No RFAM records parsed from Rfam.cm"
        )

    if clan_file is not None:

        clan_df = parse_rfam_clans(clan_file)

        if not clan_df.empty:

            rfam_df = rfam_df.merge(
                clan_df,
                on="family",
                how="left"
            )

    return rfam_df


# ============================================================
# Lookup dictionaries
# ============================================================

def build_family_name_dict(rfam_df):

    return dict(zip(rfam_df["family"], rfam_df["name"]))


def build_family_clen_dict(rfam_df):

    return dict(zip(rfam_df["family"], rfam_df["clen"]))


def build_family_nseq_dict(rfam_df):

    return dict(zip(rfam_df["family"], rfam_df["nseq"]))


def build_family_clan_dict(rfam_df):

    if "clan" not in rfam_df.columns:
        return {}

    return dict(zip(rfam_df["family"], rfam_df["clan"]))


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--cm",
        required=True,
        help="Path to Rfam.cm"
    )

    parser.add_argument(
        "--clans",
        default=None,
        help="Path to Rfam.clanin"
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV"
    )

    args = parser.parse_args()

    df = build_rfam_metadata(
        args.cm,
        args.clans
    )

    df.to_csv(
        args.out,
        sep="\t",
        index=False
    )

    print(f"\nSaved metadata table -> {args.out}")
