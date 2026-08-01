# plot_phase1_figures.py

"""
Full Phase 1 plotting pipeline for the RNASSTR revision.

Generates ALL figures that do NOT require:
- model F1 scores
- MCC scores
- inference outputs

Includes:
- KDE plots
- linear/log-axis variants
- CLEN analyses
- NSEQ analyses
- UMAPs
- split summaries

Uses:
- sequence_metadata.tsv
- family_metadata.tsv

Outputs:
- SVG + PNG for every figure

Style:
- Feng Zhang-inspired palette
- Arial
- reviewer-friendly smoothing
"""

from pathlib import Path
import argparse

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde

from sklearn.preprocessing import StandardScaler

import umap

from shared.style import (
    SPLIT_COLORS,
    setup_axis,
    save_figure,
    KDE_ALPHA,
    LINEWIDTH,
    DEFAULT_FIGSIZE,
    FENG_PINK,
)

# ============================================================
# KDE bandwidth
# ============================================================

KDE_BW = 0.25

# ============================================================
# KDE helper
# ============================================================

def kde_plot(
    ax,
    values,
    color,
    label,
    log_x=False,
    normalize=False
):

    values = np.asarray(values)

    values = values[np.isfinite(values)]

    values = values[values > 0]

    if len(values) < 10:
        return

    kde = gaussian_kde(
        values,
        bw_method=KDE_BW
    )

    xs = np.linspace(
        values.min(),
        values.max(),
        500
    )

    ys = kde(xs)

    if normalize:
        ys = ys / ys.max()

    ax.plot(
        xs,
        ys,
        color=color,
        linewidth=LINEWIDTH,
        label=label
    )

    ax.fill_between(
        xs,
        ys,
        color=color,
        alpha=KDE_ALPHA
    )

    if log_x:
        ax.set_xscale("log")


# ============================================================
# Split KDE wrapper
# ============================================================

def plot_split_kde(
    df,
    column,
    out_prefix,
    xlabel,
    log_x=False,
    normalize=False,
    subsample=None
):

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    for split in ["train", "val", "test"]:

        sub = df[df["split"] == split]

        if subsample is not None:

            if len(sub) > subsample:

                sub = sub.sample(
                    subsample,
                    random_state=1
                )

        kde_plot(
            ax=ax,
            values=sub[column],
            color=SPLIT_COLORS[split],
            label=split,
            log_x=log_x,
            normalize=normalize
        )

    setup_axis(ax)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")

    ax.legend()

    save_figure(fig, out_prefix)

    plt.close()


# ============================================================
# Scatter plot
# ============================================================

def scatter_plot(
    x,
    y,
    xlabel,
    ylabel,
    out_prefix,
    color=FENG_PINK,
    logx=False
):

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    ax.scatter(
        x,
        y,
        s=10,
        alpha=0.5,
        color=color
    )

    if logx:
        ax.set_xscale("log")

    setup_axis(ax)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    save_figure(fig, out_prefix)

    plt.close()


# ============================================================
# UMAP
# ============================================================

def run_umap(
    df,
    feature_cols,
    color_col,
    out_prefix,
    title=None,
    family_name_map=None
):

    X = df[feature_cols].copy()

    X = X.replace([np.inf, -np.inf], np.nan)

    X = X.dropna()

    keep_idx = X.index

    subdf = df.loc[keep_idx]

    scaler = StandardScaler()

    X_scaled = scaler.fit_transform(X)

    reducer = umap.UMAP(
        n_neighbors=25,
        min_dist=0.1,
        metric="euclidean",
        random_state=1
    )

    embedding = reducer.fit_transform(X_scaled)

    fig, ax = plt.subplots(
        figsize=(7.8, 6) if color_col != "split" else (6, 6)
    )

    if color_col == "split":

        for split in ["train", "val", "test"]:

            idx = subdf[color_col] == split

            ax.scatter(
                embedding[idx, 0],
                embedding[idx, 1],
                s=4,
                alpha=0.4,
                color=SPLIT_COLORS[split],
                label=split
            )

    else:

        families = (
            subdf[color_col]
            .value_counts()
            .head(10)
            .index
        )

        for fam in families:

            idx = subdf[color_col] == fam

            readable_names = {
                "RF00005": "tRNA",
                "RF00001": "5S rRNA",
                "RF00017": "Metazoan SRP RNA",
                "RF00097": "Plant snoRNA R71",
                "RF00163": "Type I hammerhead ribozyme",
                "RF00026": "U6 snRNA",
                "RF00436": "UnaL2 LINE 3′ element",
                "RF00230": "T-box leader RNA",
                "RF02543": "Eukaryotic LSU rRNA",
                "RF00174": "Cobalamin riboswitch aptamer",
            }
            if fam in readable_names:
                legend_label = f"{readable_names[fam]} ({fam})"
            elif family_name_map and fam in family_name_map:
                legend_label = f"{str(family_name_map[fam]).replace('_', ' ')} ({fam})"
            else:
                legend_label = fam

            ax.scatter(
                embedding[idx, 0],
                embedding[idx, 1],
                s=5,
                alpha=0.5,
                label=legend_label
            )

    setup_axis(ax)

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")

    if title:
        ax.set_title(title)

    if color_col == "split":
        ax.legend(markerscale=3, fontsize=8)
    else:
        ax.legend(
            markerscale=3,
            fontsize=8,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            borderaxespad=0
        )
        fig.subplots_adjust(right=0.66)

    save_figure(fig, out_prefix)

    plt.close()


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sequence-metadata",
        required=True
    )

    parser.add_argument(
        "--family-metadata",
        required=True
    )

    parser.add_argument(
        "--outdir",
        required=True
    )

    parser.add_argument(
        "--family-umap-only",
        action="store_true",
        help="Regenerate only the top-family UMAP."
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Load
    # --------------------------------------------------------

    seq_df = pd.read_csv(
        args.sequence_metadata,
        sep="\t"
    )

    fam_df = pd.read_csv(
        args.family_metadata,
        sep="\t"
    )

    print("\nLoaded sequence metadata:")
    print(len(seq_df))

    print("\nLoaded family metadata:")
    print(len(fam_df))

    family_name_map = dict(
        zip(fam_df["family"], fam_df["name"])
    )

    if args.family_umap_only:
        feature_cols = [
            "length",
            "gc_content",
            "fraction_paired",
            "stem_count",
            "basepair_density",
            "num_basepairs"
        ]
        run_umap(
            seq_df.sample(
                min(50000, len(seq_df)),
                random_state=1
            ),
            feature_cols,
            "family",
            outdir / "umap_top_families_human_readable",
            title="Structural feature UMAP by family",
            family_name_map=family_name_map
        )
        print(f"\nUpdated family UMAP written to: {outdir}")
        return

    # ========================================================
    # Length KDEs
    # ========================================================

    plot_split_kde(
        seq_df,
        "length",
        outdir / "length_kde_linear",
        "Sequence length"
    )

    plot_split_kde(
        seq_df,
        "length",
        outdir / "length_kde_logx",
        "Sequence length",
        log_x=True
    )

    # ========================================================
    # GC KDEs
    # ========================================================

    plot_split_kde(
        seq_df,
        "gc_content",
        outdir / "gc_content_kde",
        "GC content"
    )

    plot_split_kde(
        seq_df,
        "gc_content",
        outdir / "gc_content_kde_subsampled_logx",
        "GC content",
        log_x=True,
        subsample=2000
    )

    # ========================================================
    # Fraction paired
    # ========================================================

    plot_split_kde(
        seq_df,
        "fraction_paired",
        outdir / "fraction_paired_kde",
        "Fraction paired"
    )

    plot_split_kde(
        seq_df,
        "fraction_paired",
        outdir / "fraction_paired_kde_subsampled_logx",
        "Fraction paired",
        log_x=True,
        subsample=2000
    )

    # ========================================================
    # Stem counts
    # ========================================================

    plot_split_kde(
        seq_df,
        "stem_count",
        outdir / "stem_count_kde_linear",
        "Stem count"
    )

    plot_split_kde(
        seq_df,
        "stem_count",
        outdir / "stem_count_kde_logx",
        "Stem count",
        log_x=True
    )

    # ========================================================
    # CLEN
    # ========================================================

    plot_split_kde(
        fam_df,
        "clen",
        outdir / "clen_kde_logx",
        "Consensus model length (CLEN)",
        log_x=True
    )

    # ========================================================
    # NSEQ
    # ========================================================

    plot_split_kde(
        fam_df,
        "nseq",
        outdir / "nseq_kde_logx",
        "Seed alignment size (NSEQ)",
        log_x=True
    )

    # ========================================================
    # Scatter plots
    # ========================================================

    scatter_plot(
        fam_df["clen"],
        fam_df["length_mean"],
        "CLEN",
        "Mean family length",
        outdir / "clen_vs_family_length"
    )

    scatter_plot(
        fam_df["nseq"],
        fam_df["gc_content_mean"],
        "NSEQ",
        "Mean GC content",
        outdir / "nseq_vs_gc",
        logx=True
    )

    scatter_plot(
        fam_df["nseq"],
        fam_df["length_mean"],
        "NSEQ",
        "Mean family length",
        outdir / "nseq_vs_length",
        logx=True
    )

    # ========================================================
    # UMAPs
    # ========================================================

    feature_cols = [

        "length",
        "gc_content",
        "fraction_paired",
        "stem_count",
        "basepair_density",
        "num_basepairs"
    ]

    run_umap(
        seq_df.sample(
            min(50000, len(seq_df)),
            random_state=1
        ),
        feature_cols,
        "split",
        outdir / "umap_by_split",
        title="Structural feature UMAP by split"
    )

    run_umap(
        seq_df.sample(
            min(50000, len(seq_df)),
            random_state=1
        ),
        feature_cols,
        "family",
        outdir / "umap_top_families",
        title="Structural feature UMAP by family",
        family_name_map=family_name_map
    )

    print("\n===================================")
    print("DONE")
    print("===================================")

    print(f"\nFigures written to:")
    print(outdir)


if __name__ == "__main__":
    main()
