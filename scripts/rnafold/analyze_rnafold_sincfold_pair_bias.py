#!/usr/bin/env python3
"""Matched RNAfold vs retrained-SincFold base-pair bias analysis.

Inputs are the finalized per-sequence scoring files. Outputs:
  - matched per-sequence metrics
  - global paired summary
  - false-positive location/category summary
  - publication-ready 2x3 SI figure (SVG and PNG)
"""

from __future__ import annotations

import argparse
import ast
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


FENG_PINK = "#FD4078"
FENG_BLUE = "#51ADE0"
INK = "#2E3440"
GRID = "#D9DEE5"


def parse_pairs(value) -> set[tuple[int, int]]:
    if pd.isna(value) or value == "":
        return set()
    pairs = ast.literal_eval(value) if isinstance(value, str) else value
    return {tuple(sorted((int(a), int(b)))) for a, b in pairs}


def bootstrap_ci(values, rng, iterations=2000):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan
    means = np.empty(iterations)
    for i in range(iterations):
        means[i] = rng.choice(x, size=len(x), replace=True).mean()
    return tuple(np.quantile(means, [0.025, 0.975]))


def safe_wilcoxon(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = np.asarray(a)[mask], np.asarray(b)[mask]
    if len(a) == 0 or np.allclose(a, b):
        return np.nan, np.nan, len(a)
    stat, p = wilcoxon(a, b, zero_method="wilcox", alternative="two-sided")
    return float(stat), float(p), len(a)


def fp_category(fp, gt_pairs, gt_partner):
    i, j = fp
    if i in gt_partner or j in gt_partner:
        return "Alternative partner"
    for a, b in gt_pairs:
        if abs(i - a) <= 2 and abs(j - b) <= 2:
            return "Adjacent to annotated pair"
    return "Other unannotated region"


def load_matched(rnafold_path: Path, sincfold_path: Path) -> pd.DataFrame:
    rcols = [
        "id", "family", "sequence", "length", "ground_truth_base_pairs",
        "prediction_base_pairs", "tp", "fp", "fn", "precision", "recall",
        "f1", "mcc", "status",
    ]
    rna = pd.read_csv(rnafold_path, usecols=rcols)
    rna = rna[rna["status"].eq("scored")].drop(columns="status")
    rna = rna.rename(columns={
        "prediction_base_pairs": "rnafold_prediction_base_pairs",
        "tp": "rnafold_tp", "fp": "rnafold_fp", "fn": "rnafold_fn",
        "precision": "rnafold_precision", "recall": "rnafold_recall",
        "f1": "rnafold_f1", "mcc": "rnafold_mcc",
    })
    wanted = set(rna["id"])

    scols = [
        "id", "retrained_status", "retrained_prediction_base_pairs",
        "retrained_tp", "retrained_fp", "retrained_fn", "retrained_precision",
        "retrained_recall", "retrained_f1", "retrained_mcc",
    ]
    chunks = []
    for chunk in pd.read_csv(sincfold_path, usecols=scols, chunksize=250_000):
        hit = chunk[chunk["id"].isin(wanted) & chunk["retrained_status"].eq("scored")]
        if not hit.empty:
            chunks.append(hit.drop(columns="retrained_status"))
    sinc = pd.concat(chunks, ignore_index=True)
    matched = rna.merge(sinc, on="id", how="inner", validate="one_to_one")
    return matched


def calculate_pair_details(matched: pd.DataFrame):
    rows = []
    fp_records = []
    for rec in matched.itertuples(index=False):
        gt = parse_pairs(rec.ground_truth_base_pairs)
        rp = parse_pairs(rec.rnafold_prediction_base_pairs)
        sp = parse_pairs(rec.retrained_prediction_base_pairs)
        gt_n = len(gt)
        length = int(rec.length)
        gt_partner = {}
        for i, j in gt:
            gt_partner[i] = j
            gt_partner[j] = i

        row = {
            "id": rec.id,
            "family": rec.family,
            "length": length,
            "ground_truth_pairs": gt_n,
        }
        for label, pairs in (("rnafold", rp), ("retrained_sincfold", sp)):
            fp_pairs = pairs - gt
            fn_pairs = gt - pairs
            pred_n = len(pairs)
            row[f"{label}_predicted_pairs"] = pred_n
            row[f"{label}_pair_difference"] = pred_n - gt_n
            row[f"{label}_pair_ratio"] = pred_n / gt_n if gt_n else np.nan
            row[f"{label}_fp_per_reference_pair"] = len(fp_pairs) / gt_n if gt_n else np.nan
            row[f"{label}_fn_per_reference_pair"] = len(fn_pairs) / gt_n if gt_n else np.nan
            row[f"{label}_precision"] = getattr(rec, f"{'rnafold' if label == 'rnafold' else 'retrained'}_precision")
            row[f"{label}_recall"] = getattr(rec, f"{'rnafold' if label == 'rnafold' else 'retrained'}_recall")
            row[f"{label}_f1"] = getattr(rec, f"{'rnafold' if label == 'rnafold' else 'retrained'}_f1")
            row[f"{label}_mcc"] = getattr(rec, f"{'rnafold' if label == 'rnafold' else 'retrained'}_mcc")
            for pair in fp_pairs:
                fp_records.append({
                    "id": rec.id,
                    "family": rec.family,
                    "model": "RNAfold" if label == "rnafold" else "Retrained SincFold",
                    "i": pair[0],
                    "j": pair[1],
                    "normalized_midpoint": ((pair[0] + pair[1]) / 2) / max(length - 1, 1),
                    "normalized_span": (pair[1] - pair[0]) / max(length - 1, 1),
                    "category": fp_category(pair, gt, gt_partner),
                })
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(fp_records)


def make_summaries(per_sequence, fp_records, outdir):
    rng = np.random.default_rng(20260724)
    metric_pairs = [
        ("Predicted base pairs", "rnafold_predicted_pairs", "retrained_sincfold_predicted_pairs"),
        ("Predicted/reference pair ratio", "rnafold_pair_ratio", "retrained_sincfold_pair_ratio"),
        ("False positives/reference pair", "rnafold_fp_per_reference_pair", "retrained_sincfold_fp_per_reference_pair"),
        ("False negatives/reference pair", "rnafold_fn_per_reference_pair", "retrained_sincfold_fn_per_reference_pair"),
        ("Precision", "rnafold_precision", "retrained_sincfold_precision"),
        ("Recall", "rnafold_recall", "retrained_sincfold_recall"),
        ("F1", "rnafold_f1", "retrained_sincfold_f1"),
        ("MCC", "rnafold_mcc", "retrained_sincfold_mcc"),
    ]
    rows = []
    for name, rcol, scol in metric_pairs:
        a = per_sequence[rcol].to_numpy(float)
        b = per_sequence[scol].to_numpy(float)
        stat, p, n = safe_wilcoxon(a, b)
        for model, x in (("RNAfold", a), ("Retrained SincFold", b)):
            x = x[np.isfinite(x)]
            lo, hi = bootstrap_ci(x, rng)
            rows.append({
                "metric": name,
                "model": model,
                "n": len(x),
                "mean": np.mean(x),
                "median": np.median(x),
                "sd": np.std(x, ddof=1),
                "mean_ci95_low": lo,
                "mean_ci95_high": hi,
                "paired_wilcoxon_n": n,
                "paired_wilcoxon_statistic": stat,
                "paired_wilcoxon_p": p,
            })
    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "rnafold_vs_retrained_pair_bias_summary.csv", index=False)

    if fp_records.empty:
        fp_summary = pd.DataFrame()
    else:
        fp_summary = (
            fp_records.groupby(["model", "category"], observed=True)
            .size().rename("false_positive_pairs").reset_index()
        )
        totals = fp_summary.groupby("model")["false_positive_pairs"].transform("sum")
        fp_summary["fraction_of_model_false_positives"] = fp_summary["false_positive_pairs"] / totals
    fp_summary.to_csv(outdir / "rnafold_vs_retrained_fp_location_summary.csv", index=False)
    return summary, fp_summary


def plot_figure(per_sequence, fp_records, summary, outdir):
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
    })
    fig, axs = plt.subplots(2, 3, figsize=(12.2, 7.4))
    axa, axb, axc, axd, axe, axf = axs.flat
    cmap_blue = LinearSegmentedColormap.from_list("fengblue", ["#EFF8FD", FENG_BLUE, "#176B98"])
    cmap_pink = LinearSegmentedColormap.from_list("fengpink", ["#FFF0F5", FENG_PINK, "#A71947"])

    maxima = np.nanpercentile(np.r_[
        per_sequence["ground_truth_pairs"],
        per_sequence["rnafold_predicted_pairs"],
        per_sequence["retrained_sincfold_predicted_pairs"],
    ], 99.5)
    maxima = max(10, math.ceil(maxima / 10) * 10)
    for ax, ycol, title, cmap in (
        (axa, "rnafold_predicted_pairs", "RNAfold", cmap_blue),
        (axb, "retrained_sincfold_predicted_pairs", "Retrained SincFold", cmap_pink),
    ):
        ax.hexbin(
            per_sequence["ground_truth_pairs"], per_sequence[ycol],
            gridsize=45, mincnt=1, cmap=cmap, bins="log", linewidths=0,
            extent=(0, maxima, 0, maxima),
        )
        ax.plot([0, maxima], [0, maxima], color=INK, lw=1.2, ls="--")
        ax.set_xlim(0, maxima)
        ax.set_ylim(0, maxima)
        ax.set_xlabel("Reference base pairs")
        ax.set_ylabel("Predicted base pairs")
        ax.set_title(title, fontsize=12, weight="bold")

    ratio_data = [
        per_sequence["rnafold_pair_ratio"].dropna().clip(upper=3),
        per_sequence["retrained_sincfold_pair_ratio"].dropna().clip(upper=3),
    ]
    parts = axc.violinplot(ratio_data, positions=[1, 2], showmedians=True, showextrema=False)
    for body, color in zip(parts["bodies"], [FENG_BLUE, FENG_PINK]):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.30)
    parts["cmedians"].set_color([FENG_BLUE, FENG_PINK])
    parts["cmedians"].set_linewidth(2)
    axc.axhline(1, color=INK, lw=1.2, ls="--")
    axc.set_xticks([1, 2], ["RNAfold", "Retrained\nSincFold"])
    axc.set_ylabel("Predicted/reference pair ratio")
    axc.set_ylim(0, 3)

    measures = ["False positives/reference pair", "False negatives/reference pair"]
    x = np.arange(2)
    width = 0.34
    for offset, model, color in ((-width/2, "RNAfold", FENG_BLUE), (width/2, "Retrained SincFold", FENG_PINK)):
        vals, lo, hi = [], [], []
        for metric in measures:
            row = summary[(summary.metric == metric) & (summary.model == model)].iloc[0]
            vals.append(row["mean"])
            lo.append(row["mean"] - row["mean_ci95_low"])
            hi.append(row["mean_ci95_high"] - row["mean"])
        axd.bar(x + offset, vals, width, facecolor=mpl.colors.to_rgba(color, 0.22), edgecolor=color, linewidth=1.8, label=model)
        axd.errorbar(x + offset, vals, yerr=[lo, hi], fmt="none", ecolor=color, elinewidth=1.2, capsize=3)
    axd.set_xticks(x, ["False positive", "False negative"])
    axd.set_ylabel("Errors per reference pair")
    axd.legend(loc="upper right")

    cats = ["Alternative partner", "Adjacent to annotated pair", "Other unannotated region"]
    fp_counts = fp_records.groupby(["model", "category"]).size().unstack(fill_value=0).reindex(
        index=["RNAfold", "Retrained SincFold"], columns=cats, fill_value=0
    )
    fp_frac = fp_counts.div(fp_counts.sum(axis=1), axis=0)
    bottoms = np.zeros(2)
    cat_colors = [INK, "#9FA8B3", "#D9DEE5"]
    for cat, color in zip(cats, cat_colors):
        vals = fp_frac[cat].to_numpy()
        axe.bar([0, 1], vals, bottom=bottoms, color=color, width=0.62, label=cat)
        bottoms += vals
    axe.set_xticks([0, 1], ["RNAfold", "Retrained\nSincFold"])
    axe.set_ylabel("Fraction of false-positive pairs")
    axe.set_ylim(0, 1)
    axe.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.02))

    bins = np.linspace(0, 1, 21)
    centers = (bins[:-1] + bins[1:]) / 2
    for model, color in (("RNAfold", FENG_BLUE), ("Retrained SincFold", FENG_PINK)):
        vals = fp_records.loc[fp_records.model.eq(model), "normalized_midpoint"].to_numpy()
        hist, _ = np.histogram(vals, bins=bins, density=True)
        axf.plot(centers, hist, color=color, lw=2, label=model)
    axf.set_xlabel("Normalized sequence position")
    axf.set_ylabel("False-positive pair density")
    axf.set_xlim(0, 1)
    axf.legend()

    for label, ax in zip("ABCDEF", axs.flat):
        ax.text(-0.18, 1.08, label, transform=ax.transAxes, fontsize=17, weight="bold", va="top")
        ax.grid(axis="y", color=GRID, lw=0.6, alpha=0.7)
        ax.tick_params(labelsize=9)

    fig.tight_layout(w_pad=2.0, h_pad=2.2)
    prefix = outdir / "supplementary_pair_bias_rnafold_vs_retrained_sincfold"
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)

    # Compact version containing only the base-pair count comparisons.
    fig_ab, (axa, axb) = plt.subplots(1, 2, figsize=(8.2, 4.0), sharex=True, sharey=True)
    for ax, ycol, title, cmap in (
        (axa, "rnafold_predicted_pairs", "RNAfold", cmap_blue),
        (axb, "retrained_sincfold_predicted_pairs", "Retrained SincFold", cmap_pink),
    ):
        ax.hexbin(
            per_sequence["ground_truth_pairs"], per_sequence[ycol],
            gridsize=45, mincnt=1, cmap=cmap, bins="log", linewidths=0,
            extent=(0, maxima, 0, maxima),
        )
        ax.plot([0, maxima], [0, maxima], color=INK, lw=1.2, ls="--")
        ax.set_xlim(0, maxima)
        ax.set_ylim(0, maxima)
        ax.set_xlabel("Reference base pairs")
        ax.set_ylabel("Predicted base pairs")
        ax.set_title(title, fontsize=12, weight="bold")
        ax.grid(axis="y", color=GRID, lw=0.6, alpha=0.7)
        ax.tick_params(labelsize=9)
    for label, ax in zip("AB", (axa, axb)):
        ax.text(-0.17, 1.08, label, transform=ax.transAxes, fontsize=17, weight="bold", va="top")
    fig_ab.tight_layout(w_pad=2.0)
    prefix_ab = outdir / "supplementary_pair_counts_rnafold_vs_retrained_sincfold_AB"
    fig_ab.savefig(prefix_ab.with_suffix(".svg"), bbox_inches="tight")
    fig_ab.savefig(prefix_ab.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig_ab)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rnafold", type=Path, required=True)
    parser.add_argument("--sincfold", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    matched = load_matched(args.rnafold, args.sincfold)
    per_sequence, fp_records = calculate_pair_details(matched)
    per_sequence.to_csv(args.outdir / "rnafold_vs_retrained_pair_bias_per_sequence.csv", index=False)
    fp_records.to_csv(args.outdir / "rnafold_vs_retrained_false_positive_pairs.csv", index=False)
    summary, _ = make_summaries(per_sequence, fp_records, args.outdir)
    plot_figure(per_sequence, fp_records, summary, args.outdir)

    report = [
        f"Matched scored sequences: {len(per_sequence):,}",
        f"Families: {per_sequence['family'].nunique():,}",
        f"Sequences with >=1 reference pair: {(per_sequence['ground_truth_pairs'] > 0).sum():,}",
        f"RNAfold false-positive pairs: {(fp_records['model'] == 'RNAfold').sum():,}",
        f"Retrained SincFold false-positive pairs: {(fp_records['model'] == 'Retrained SincFold').sum():,}",
    ]
    (args.outdir / "rnafold_vs_retrained_pair_bias_report.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    main()
