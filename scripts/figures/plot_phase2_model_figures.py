#!/usr/bin/env python3
"""Generate Phase 2 SI figures from standardized family-level model summaries.

Models included:
  - published SincFold
  - retrained SincFold
  - Lyra TransPred

RNAfold is intentionally excluded from this family-level comparison pipeline.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess


FENG_PINK = "#FD4078"
FENG_GREY = "#6E6E6E"
FENG_BLUE = "#51ADE0"
FENG_GOLD = "#E6A53A"

MODEL_STYLE = {
    "published_sincfold": {
        "label": "Published SincFold",
        "color": FENG_GREY,
        "marker": "o",
    },
    "retrained_sincfold": {
        "label": "Retrained SincFold",
        "color": FENG_PINK,
        "marker": "s",
    },
    "lyra_transpred": {
        "label": "Lyra-TransPred",
        "color": FENG_BLUE,
        "marker": "^",
    },
}

mpl.rcParams.update({
    "font.family": "Arial",
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def setup_axis(ax):
    ax.grid(True, alpha=0.2, linewidth=0.8)
    ax.tick_params(axis="both", labelsize=10)


def save_figure(fig, prefix: Path):
    fig.savefig(prefix.with_suffix(".png"))
    fig.savefig(prefix.with_suffix(".svg"))
    plt.close(fig)


def load_sincfold(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {
        "family", "model", "n_truth", "n_scored", "coverage",
        "mean_f1", "mean_mcc",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    df = df[df["model"].isin(["published", "retrained"])].copy()
    df["model"] = df["model"].map({
        "published": "published_sincfold",
        "retrained": "retrained_sincfold",
    })
    return df


def load_transpred(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {
        "family", "model", "n_truth", "n_scored", "coverage",
        "mean_f1", "mean_mcc",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    df = df[df["model"] == "transpred"].copy()
    df["model"] = "lyra_transpred"
    return df


def load_family_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"family", "split", "stem_count_mean"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    df = df[df["split"] == "test"].copy()
    name_column = "name" if "name" in df.columns else None
    if name_column:
        df["family_label"] = np.where(
            df[name_column].notna() & (df[name_column].astype(str) != ""),
            df[name_column].astype(str) + " (" + df["family"].astype(str) + ")",
            df["family"].astype(str),
        )
    else:
        df["family_label"] = df["family"]
    return df[["family", "family_label", "stem_count_mean"]]


def pair_models(
    family_scores: pd.DataFrame, model_x: str, model_y: str
) -> pd.DataFrame:
    columns = [
        "family", "n_truth", "n_scored", "coverage", "mean_f1", "mean_mcc",
    ]
    x = family_scores[family_scores["model"] == model_x][columns].copy()
    y = family_scores[family_scores["model"] == model_y][columns].copy()
    merged = x.merge(y, on="family", suffixes=("_x", "_y"), how="inner")
    return merged.dropna(subset=[
        "mean_f1_x", "mean_f1_y", "mean_mcc_x", "mean_mcc_y",
    ])


def point_sizes(n: pd.Series) -> np.ndarray:
    values = np.maximum(pd.to_numeric(n, errors="coerce").fillna(1), 1)
    return 14 + 8 * np.log10(values)


def comparison_panel(
    ax,
    data: pd.DataFrame,
    metric: str,
    x_label: str,
    y_label: str,
    color: str,
    panel_label: str,
):
    x = data[f"{metric}_x"]
    y = data[f"{metric}_y"]
    sizes = point_sizes(np.minimum(data["n_scored_x"], data["n_scored_y"]))
    ax.scatter(
        x, y, s=sizes, color=color, alpha=0.55,
        edgecolors="none", rasterized=True,
    )
    ax.plot([0, 1], [0, 1], linestyle="--", color=FENG_GREY, linewidth=1.2)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    setup_axis(ax)
    rho, pvalue = spearmanr(x, y)
    delta = np.median(y - x)
    ax.text(
        -0.12, 1.05, panel_label,
        transform=ax.transAxes, fontsize=13, fontweight="bold",
    )
    return {
        "n_families": len(data),
        "spearman_rho": rho,
        "spearman_pvalue": pvalue,
        "median_delta_y_minus_x": delta,
    }


def make_comparison_figure(
    data: pd.DataFrame,
    model_x: str,
    model_y: str,
    output_prefix: Path,
    point_color: str,
):
    style_x = MODEL_STYLE[model_x]
    style_y = MODEL_STYLE[model_y]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.1))
    stats = []
    stats.append({
        "metric": "F1",
        **comparison_panel(
            axes[0], data, "mean_f1",
            f"{style_x['label']} mean family F1",
            f"{style_y['label']} mean family F1",
            point_color, "A",
        ),
    })
    stats.append({
        "metric": "MCC",
        **comparison_panel(
            axes[1], data, "mean_mcc",
            f"{style_x['label']} mean family MCC",
            f"{style_y['label']} mean family MCC",
            point_color, "B",
        ),
    })
    fig.tight_layout(w_pad=2.5)
    save_figure(fig, output_prefix)
    return pd.DataFrame(stats)


def complexity_panel(ax, data, metric, panel_label):
    for model in [
        "published_sincfold", "retrained_sincfold", "lyra_transpred"
    ]:
        sub = data[
            (data["model"] == model)
            & data[metric].notna()
            & (data["stem_count_mean"] > 0)
        ].copy()
        style = MODEL_STYLE[model]
        sizes = point_sizes(sub["n_scored"])
        ax.scatter(
            sub["stem_count_mean"], sub[metric],
            s=sizes, alpha=0.24, color=style["color"],
            marker=style["marker"], edgecolors="none", rasterized=True,
            label=style["label"],
        )
        trend_input = (
            sub[["stem_count_mean", metric]]
            .dropna()
            .sort_values("stem_count_mean")
        )
        if len(trend_input) >= 20:
            trend = lowess(
                trend_input[metric],
                trend_input["stem_count_mean"],
                frac=0.35,
                it=1,
                return_sorted=True,
            )
            trend[:, 1] = np.clip(trend[:, 1], 0, 1)
            ax.plot(
                trend[:, 0], trend[:, 1],
                color=style["color"], linewidth=2.2,
            )
    ax.set_xscale("log")
    ax.set_xlim(0.8, max(60, data["stem_count_mean"].max() * 1.08))
    ax.set_xticks([1, 2, 5, 10, 20, 50])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Mean annotated stem count per family")
    ax.set_ylabel("Mean family F1" if metric == "mean_f1" else "Mean family MCC")
    setup_axis(ax)
    ax.text(
        -0.12, 1.05, panel_label,
        transform=ax.transAxes, fontsize=13, fontweight="bold",
    )



def make_complexity_figure(data: pd.DataFrame, output_prefix: Path):
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0), sharex=True)
    complexity_panel(axes[0], data, "mean_f1", "A")
    complexity_panel(axes[1], data, "mean_mcc", "B")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.04),
        ncol=3, fontsize=9,
    )
    fig.tight_layout(w_pad=2.2)
    save_figure(fig, output_prefix)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sincfold-family", required=True, type=Path)
    parser.add_argument("--lyra-family", required=True, type=Path)
    parser.add_argument("--family-metadata", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    sincfold = load_sincfold(args.sincfold_family)
    transpred = load_transpred(args.lyra_family)
    metadata = load_family_metadata(args.family_metadata)
    all_scores = pd.concat([sincfold, transpred], ignore_index=True)

    sincfold_pair = pair_models(
        all_scores, "published_sincfold", "retrained_sincfold"
    )
    lyra_pair = pair_models(
        all_scores, "published_sincfold", "lyra_transpred"
    )

    sincfold_stats = make_comparison_figure(
        sincfold_pair,
        "published_sincfold", "retrained_sincfold",
        args.outdir / "si_sincfold_family_comparison",
        FENG_PINK,
    )
    lyra_stats = make_comparison_figure(
        lyra_pair,
        "published_sincfold", "lyra_transpred",
        args.outdir / "si_transpred_vs_published_sincfold",
        FENG_BLUE,
    )

    complexity = all_scores.merge(metadata, on="family", how="inner")
    make_complexity_figure(
        complexity, args.outdir / "si_performance_vs_stem_count"
    )

    sincfold_pair.to_csv(
        args.outdir / "sincfold_family_comparison.csv", index=False
    )
    lyra_pair.to_csv(
        args.outdir / "transpred_vs_published_family_comparison.csv", index=False
    )
    complexity.to_csv(
        args.outdir / "complexity_model_metrics.csv", index=False
    )
    top_complex_families = (
        metadata.nlargest(10, "stem_count_mean")["family"].tolist()
    )
    (
        complexity[complexity["family"].isin(top_complex_families)]
        .sort_values(["stem_count_mean", "family", "model"], ascending=[False, True, True])
        .to_csv(args.outdir / "top_complex_family_metrics.csv", index=False)
    )
    stats = pd.concat([
        sincfold_stats.assign(comparison="published_vs_retrained_sincfold"),
        lyra_stats.assign(comparison="published_sincfold_vs_lyra_transpred"),
    ], ignore_index=True)
    stats.to_csv(args.outdir / "phase2_figure_statistics.csv", index=False)

    coverage = (
        all_scores.groupby("model", as_index=False)
        .agg(
            families=("family", "nunique"),
            truth_sequences=("n_truth", "sum"),
            scored_sequences=("n_scored", "sum"),
        )
    )
    coverage["coverage"] = (
        coverage["scored_sequences"] / coverage["truth_sequences"]
    )
    coverage.to_csv(args.outdir / "phase2_coverage_summary.csv", index=False)

    print(f"Figures and tables written to {args.outdir}")


if __name__ == "__main__":
    main()
