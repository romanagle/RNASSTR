#!/usr/bin/env python3
"""Rebuild Supplementary Figures 1–4 in the RNASSTR Feng-style palette."""

from __future__ import annotations

import argparse
import ast
import csv
from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


FENG_PINK = "#FD4078"
FENG_GREY = "#6E6E6E"
FENG_BLUE = "#51ADE0"
FENG_GOLD = "#E6A53A"

mpl.rcParams.update({
    "font.family": "Arial",
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

FENG_BLUE_CMAP = LinearSegmentedColormap.from_list(
    "feng_blue", ["#FFFFFF", "#D7EEF9", FENG_BLUE, "#146B98"]
)


def setup_axis(ax):
    ax.grid(True, alpha=0.2, linewidth=0.8)
    ax.tick_params(axis="both", labelsize=10)


def panel_label(ax, label):
    ax.text(
        -0.13, 1.05, label, transform=ax.transAxes,
        fontsize=13, fontweight="bold",
    )


def save_figure(fig, prefix: Path):
    fig.savefig(prefix.with_suffix(".png"))
    fig.savefig(prefix.with_suffix(".svg"))
    plt.close(fig)


def get_length_counts(sequence_metadata: Path, cache: Path) -> pd.DataFrame:
    if cache.exists():
        return pd.read_csv(cache)
    counts = Counter()
    for chunk in pd.read_csv(
        sequence_metadata,
        sep="\t",
        usecols=["length"],
        chunksize=500_000,
    ):
        values = pd.to_numeric(chunk["length"], errors="coerce").dropna().astype(int)
        counts.update(values.tolist())
    result = pd.DataFrame(
        sorted(counts.items()), columns=["length", "count"]
    )
    result.to_csv(cache, index=False)
    return result


def make_figure1(length_counts: pd.DataFrame, outdir: Path):
    lengths = length_counts["length"].to_numpy(dtype=float)
    counts = length_counts["count"].to_numpy(dtype=float)
    positive = lengths > 0
    lengths, counts = lengths[positive], counts[positive]
    # Estimate density in log-length space, using sequence counts as weights,
    # then transform the density back to the original length scale.
    log_lengths = np.log10(lengths)
    kde = gaussian_kde(log_lengths, weights=counts, bw_method=0.075)
    grid = np.geomspace(lengths.min(), lengths.max(), 700)
    density = kde(np.log10(grid)) / (grid * np.log(10))

    fig, ax = plt.subplots(figsize=(6.3, 4.7))
    ax.plot(grid, density, color=FENG_PINK, linewidth=2.2)
    ax.fill_between(
        grid, density,
        np.full(grid.size, max(density.min(), np.finfo(float).tiny)),
        color=FENG_PINK, alpha=0.22,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Sequence length (nt)")
    ax.set_ylabel("Density (log scale)")
    setup_axis(ax)
    fig.tight_layout()
    save_figure(fig, outdir / "supplementary_figure_1_length_distribution")


def make_figure2(outdir: Path):
    # Approximate values supplied for the expanded training set. SincFold
    # values are the original measurements increased by 30%; Lyra training is
    # the midpoint of the reported 10–12 GPU h range. Its inference estimate
    # uses the same relative speed as training.
    sincfold_training = 64.0 * 1.30
    sincfold_inference = 2.0 * 1.30
    lyra_training = 11.0
    relative_speed = sincfold_training / lyra_training
    lyra_inference = sincfold_inference / relative_speed

    stages = ["Training per epoch", "Test-set inference"]
    sincfold = [sincfold_training, sincfold_inference]
    lyra = [lyra_training, lyra_inference]
    x = np.arange(len(stages))
    width = 0.34
    fig, ax = plt.subplots(figsize=(6.0, 4.7))
    bars_sincfold = ax.bar(
        x - width / 2, sincfold, width,
        color="#FFD9E5", edgecolor=FENG_PINK, linewidth=1.8,
        label="SincFold",
    )
    bars_lyra = ax.bar(
        x + width / 2, lyra, width,
        color="#DDF1FA", edgecolor=FENG_BLUE, linewidth=1.8,
        label="Lyra-TransPred",
    )
    for bars, values in [(bars_sincfold, sincfold), (bars_lyra, lyra)]:
        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value * 1.18,
                f"{value:.1f}" if value >= 1 else f"{value:.2f}",
                ha="center", va="bottom", fontsize=9,
            )
    ax.set_xticks(x, stages)
    ax.set_ylabel("GPU hours")
    ax.set_yscale("log")
    ax.set_ylim(0.2, 150)
    ax.legend(loc="upper right")
    ax.grid(False)
    ax.grid(axis="y", which="major", alpha=0.22, linewidth=0.9)
    ax.tick_params(axis="both", labelsize=10)
    fig.tight_layout()
    save_figure(fig, outdir / "supplementary_figure_2_model_timing")

    pd.DataFrame([
        {
            "model": "SincFold", "stage": "training_per_epoch",
            "gpu_hours": sincfold_training,
            "basis": "Original 64 GPU h increased by 30%",
        },
        {
            "model": "SincFold", "stage": "test_set_inference",
            "gpu_hours": sincfold_inference,
            "basis": "Original 2 GPU h increased by 30%",
        },
        {
            "model": "Lyra-TransPred", "stage": "training_per_epoch",
            "gpu_hours": lyra_training,
            "basis": "Midpoint of estimated 10–12 GPU h range",
        },
        {
            "model": "Lyra-TransPred", "stage": "test_set_inference",
            "gpu_hours": lyra_inference,
            "basis": "Estimated using the training-time speed ratio",
        },
    ]).to_csv(outdir / "supplementary_figure_2_timing_values.csv", index=False)


def get_feature_score_data(scores_path: Path, cache: Path) -> pd.DataFrame:
    if cache.exists():
        return pd.read_csv(cache)
    rows = []
    usecols = [
        "sequence", "length", "ground_truth_dot_bracket",
        "retrained_status", "retrained_f1",
    ]
    for chunk in pd.read_csv(scores_path, usecols=usecols, chunksize=50_000):
        chunk = chunk[chunk["retrained_status"] == "scored"].copy()
        structure = chunk["ground_truth_dot_bracket"].fillna("").astype(str)
        sequence = chunk["sequence"].fillna("").astype(str).str.upper()
        paired_count = structure.str.count(r"[\(\)\[\]\{\}<>]")
        chunk["num_basepairs"] = paired_count / 2
        chunk["fraction_paired"] = paired_count / chunk["length"]
        chunk["gc_content"] = (
            sequence.str.count("G") + sequence.str.count("C")
        ) / chunk["length"]
        rows.append(chunk[[
            "retrained_f1", "num_basepairs", "gc_content",
            "fraction_paired", "length",
        ]])
    result = pd.concat(rows, ignore_index=True)
    result.to_csv(cache, index=False)
    return result


def hex_panel(
    ax, x, y, xlabel, ylabel, panel,
    ytick_values=None, ytick_labels=None,
):
    clean = pd.DataFrame({"x": x, "y": y}).replace(
        [np.inf, -np.inf], np.nan
    ).dropna()
    hb = ax.hexbin(
        clean["x"], clean["y"],
        gridsize=55,
        mincnt=1,
        cmap=FENG_BLUE_CMAP,
        norm=LogNorm(),
        linewidths=0,
    )
    ax.set_xlim(-0.02, 1.02)
    if ytick_values is not None:
        ax.set_yticks(ytick_values)
        ax.set_yticklabels(ytick_labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    setup_axis(ax)
    panel_label(ax, panel)
    return hb


def make_feature_figure(features: pd.DataFrame, output_prefix: Path, xlabel: str):
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 8.2))
    pair_ticks = np.array([0, 1, 5, 10, 50, 100, 500, 1000])
    length_ticks = np.array([20, 50, 100, 200, 500, 1000, 2000, 5000])
    panels = [
        (
            np.log10(features["num_basepairs"] + 1),
            "Number of ground-truth base pairs", "A",
            np.log10(pair_ticks + 1), [f"{v:g}" for v in pair_ticks],
        ),
        (features["gc_content"], "GC content", "B", None, None),
        (features["fraction_paired"], "Fraction paired", "C", None, None),
        (
            np.log10(features["length"]),
            "Sequence length (nt)", "D",
            np.log10(length_ticks), [f"{v:g}" for v in length_ticks],
        ),
    ]
    for ax, (y, ylabel, label, tick_values, tick_labels) in zip(
        axes.flat, panels
    ):
        hb = hex_panel(
            ax, features["retrained_f1"], y,
            xlabel, ylabel, label,
            tick_values, tick_labels,
        )
        cbar = fig.colorbar(hb, ax=ax, pad=0.02)
        cbar.set_label("Sequences per bin")
    fig.tight_layout(w_pad=2.1, h_pad=2.0)
    save_figure(fig, output_prefix)


def make_figure3(features: pd.DataFrame, outdir: Path):
    make_feature_figure(
        features,
        outdir / "supplementary_figure_3_sincfold_features_vs_f1",
        "Retrained SincFold F1",
    )


def get_lyra_feature_data(scores_path: Path, cache: Path) -> pd.DataFrame:
    if cache.exists():
        return pd.read_csv(cache)
    rows = []
    usecols = [
        "sequence", "length", "ground_truth_dot_bracket",
        "transpred_status", "transpred_f1",
    ]
    for chunk in pd.read_csv(scores_path, usecols=usecols, chunksize=50_000):
        chunk = chunk[chunk["transpred_status"] == "scored"].copy()
        structure = chunk["ground_truth_dot_bracket"].fillna("").astype(str)
        sequence = chunk["sequence"].fillna("").astype(str).str.upper()
        paired_count = structure.str.count(r"[\(\)\[\]\{\}<>]")
        chunk["num_basepairs"] = paired_count / 2
        chunk["fraction_paired"] = paired_count / chunk["length"]
        chunk["gc_content"] = (
            sequence.str.count("G") + sequence.str.count("C")
        ) / chunk["length"]
        rows.append(chunk[[
            "transpred_f1", "num_basepairs", "gc_content",
            "fraction_paired", "length",
        ]])
    result = pd.concat(rows, ignore_index=True)
    result.to_csv(cache, index=False)
    return result


def make_figure4(features: pd.DataFrame, outdir: Path):
    renamed = features.rename(columns={"transpred_f1": "retrained_f1"})
    # Use the same layout, binning, scales, and color mapping as Figure 3.
    make_feature_figure(
        renamed,
        outdir / "supplementary_figure_4_lyra_transpred_features_vs_f1",
        "Lyra-TransPred F1",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-metadata", required=True, type=Path)
    parser.add_argument("--retrained-scores", required=True, type=Path)
    parser.add_argument("--lyra-scores", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    lengths = get_length_counts(
        args.sequence_metadata, args.outdir / "length_counts.csv"
    )
    make_figure1(lengths, args.outdir)
    make_figure2(args.outdir)
    features = get_feature_score_data(
        args.retrained_scores, args.outdir / "retrained_feature_scores.csv"
    )
    make_figure3(features, args.outdir)
    lyra_features = get_lyra_feature_data(
        args.lyra_scores, args.outdir / "lyra_transpred_feature_scores.csv"
    )
    make_figure4(lyra_features, args.outdir)
    print(f"Figures written to {args.outdir}")


if __name__ == "__main__":
    main()
