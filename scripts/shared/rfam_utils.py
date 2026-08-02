"""Parse Rfam covariance-model and clan metadata."""

from pathlib import Path

import pandas as pd


def parse_rfam_cm(cm_file):
    """Return ACC, NAME, CLEN, and NSEQ fields from an Rfam CM file."""
    records = []
    current = {}

    with Path(cm_file).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line == "//":
                if "family" in current:
                    records.append(current.copy())
                current = {}
                continue

            field, _, value = line.partition(" ")
            value = value.strip()
            if field == "ACC" and value:
                current["family"] = value
            elif field == "NAME" and value:
                current["name"] = value
            elif field in {"CLEN", "NSEQ"} and value:
                try:
                    current[field.lower()] = int(value)
                except ValueError:
                    continue

    dataframe = pd.DataFrame(records)
    if dataframe.empty:
        return dataframe

    if "clen" not in dataframe.columns:
        dataframe["clen"] = pd.NA
    dataframe["clen_missing"] = dataframe["clen"].isna()
    return (
        dataframe.sort_values("clen_missing")
        .drop_duplicates(subset=["family"], keep="first")
        .drop(columns=["clen_missing"])
        .reset_index(drop=True)
    )


def parse_rfam_clans(clan_file):
    """Return Rfam family-to-clan assignments from an Rfam clan file."""
    records = []
    with Path(clan_file).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 2:
                records.append({"family": fields[1], "clan": fields[0]})
    return pd.DataFrame(records, columns=["family", "clan"])


def build_rfam_metadata(cm_file, clan_file=None):
    """Combine Rfam covariance-model metadata with optional clan assignments."""
    metadata = parse_rfam_cm(cm_file)
    if metadata.empty:
        raise ValueError("No Rfam records were parsed from the CM file")
    if clan_file is not None:
        clans = parse_rfam_clans(clan_file)
        if not clans.empty:
            metadata = metadata.merge(clans, on="family", how="left")
    return metadata


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


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Build an Rfam family metadata table from CM and clan files."
    )
    parser.add_argument("--cm", required=True, help="Path to Rfam.cm")
    parser.add_argument("--clans", help="Path to Rfam.clanin")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    metadata = build_rfam_metadata(args.cm, args.clans)
    metadata.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
