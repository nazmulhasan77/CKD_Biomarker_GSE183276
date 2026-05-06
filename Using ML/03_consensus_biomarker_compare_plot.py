"""
Create consensus biomarker comparison plots.

Consensus definition:
  Genes present in all three evidence sources:
    1. Traditional CKD-up differential expression table
    2. ML CKD-up candidate table
    3. Random Forest feature-importance table

Outputs:
  Using ML/results/consensus_traditional_ml_rf_biomarkers.csv
  Using ML/figures/consensus_biomarker_score_comparison.png
  Using ML/figures/consensus_biomarker_rank_comparison.png
  Using ML/figures/consensus_biomarker_score_heatmap.png
  Using ML/figures/consensus_overlap_counts.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare consensus biomarkers from DE, ML, and Random Forest."
    )
    parser.add_argument(
        "--traditional",
        default="results/ckd_upregulated_candidate_biomarkers.csv",
        help="Traditional CKD-up DE result table.",
    )
    parser.add_argument(
        "--ml",
        default="Using ML/results/ml_upregulated_candidate_biomarkers.csv",
        help="ML CKD-up candidate table.",
    )
    parser.add_argument(
        "--rf",
        default="Using ML/results/random_forest_gene_importance.csv",
        help="Random Forest feature-importance table.",
    )
    parser.add_argument("--outdir", default="Using ML/results")
    parser.add_argument("--figdir", default="Using ML/figures")
    parser.add_argument(
        "--top-n",
        type=int,
        default=25,
        help="Number of consensus genes to show in plots.",
    )
    return parser.parse_args()


def minmax_scale(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    low = values.min()
    high = values.max()
    if not np.isfinite(low) or not np.isfinite(high) or high == low:
        return pd.Series(np.ones(len(values)), index=values.index)
    return (values - low) / (high - low)


def load_and_merge(
    traditional_path: Path, ml_path: Path, rf_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    traditional = pd.read_csv(traditional_path).copy()
    ml = pd.read_csv(ml_path).copy()
    rf = pd.read_csv(rf_path).copy()

    traditional = traditional.sort_values(
        ["rank_score", "adj_p_value"], ascending=[False, True]
    ).reset_index(drop=True)
    ml = ml.sort_values("combined_ml_biomarker_score", ascending=False).reset_index(
        drop=True
    )
    rf = rf.sort_values("random_forest_importance", ascending=False).reset_index(
        drop=True
    )

    traditional["traditional_rank"] = np.arange(1, len(traditional) + 1)
    ml["ml_rank"] = np.arange(1, len(ml) + 1)
    rf["rf_rank"] = np.arange(1, len(rf) + 1)

    traditional_cols = [
        "gene",
        "rank_score",
        "adj_p_value",
        "diff_mean_log_expr",
        "cohen_d",
        "pct_expr_case",
        "pct_expr_control",
        "traditional_rank",
    ]
    ml_cols = [
        "gene",
        "combined_ml_biomarker_score",
        "knn_single_gene_balanced_accuracy_mean",
        "knn_single_gene_roc_auc_mean",
        "ml_rank",
    ]
    rf_cols = ["gene", "random_forest_importance", "rf_rank"]

    consensus = traditional[traditional_cols].merge(ml[ml_cols], on="gene", how="inner")
    consensus = consensus.merge(rf[rf_cols], on="gene", how="inner")

    consensus["traditional_score_scaled"] = minmax_scale(consensus["rank_score"])
    consensus["ml_score_scaled"] = minmax_scale(
        consensus["combined_ml_biomarker_score"]
    )
    consensus["rf_score_scaled"] = minmax_scale(consensus["random_forest_importance"])
    consensus["consensus_score"] = consensus[
        ["traditional_score_scaled", "ml_score_scaled", "rf_score_scaled"]
    ].mean(axis=1)
    consensus["mean_rank"] = consensus[
        ["traditional_rank", "ml_rank", "rf_rank"]
    ].mean(axis=1)
    consensus = consensus.sort_values(
        ["consensus_score", "mean_rank"], ascending=[False, True]
    ).reset_index(drop=True)

    overlap_rows = []
    for n in [25, 50, 100, 200, 500]:
        n_trad = min(n, len(traditional))
        n_ml = min(n, len(ml))
        n_rf = min(n, len(rf))
        overlap = (
            set(traditional.head(n_trad)["gene"])
            & set(ml.head(n_ml)["gene"])
            & set(rf.head(n_rf)["gene"])
        )
        overlap_rows.append({"top_n": n, "triple_overlap_genes": len(overlap)})
    overlap_rows.append({"top_n": "all", "triple_overlap_genes": len(consensus)})
    overlap_counts = pd.DataFrame(overlap_rows)
    return consensus, overlap_counts


def plot_score_comparison(top: pd.DataFrame, figdir: Path) -> None:
    plot_df = top[
        [
            "gene",
            "traditional_score_scaled",
            "ml_score_scaled",
            "rf_score_scaled",
        ]
    ].copy()
    plot_df = plot_df.iloc[::-1]
    y = np.arange(len(plot_df))
    height = 0.24

    fig, ax = plt.subplots(figsize=(9, max(5.5, len(plot_df) * 0.28)), dpi=180)
    ax.barh(
        y - height,
        plot_df["traditional_score_scaled"],
        height=height,
        label="Traditional DE",
        color="#2f6f8f",
    )
    ax.barh(
        y,
        plot_df["ml_score_scaled"],
        height=height,
        label="Combined ML",
        color="#b24a3b",
    )
    ax.barh(
        y + height,
        plot_df["rf_score_scaled"],
        height=height,
        label="Random Forest",
        color="#6f8f3a",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"])
    ax.set_xlabel("Scaled score within consensus genes")
    ax.set_title("Consensus CKD-Up Biomarkers: Method Score Comparison")
    ax.legend(frameon=False, loc="lower right")
    ax.set_xlim(0, 1.05)
    fig.tight_layout()
    fig.savefig(figdir / "consensus_biomarker_score_comparison.png")
    plt.close(fig)


def plot_rank_comparison(top: pd.DataFrame, figdir: Path) -> None:
    plot_df = top.sort_values("consensus_score", ascending=True).copy()
    y = np.arange(len(plot_df))

    fig, ax = plt.subplots(figsize=(9, max(5.5, len(plot_df) * 0.28)), dpi=180)
    ax.scatter(plot_df["traditional_rank"], y, label="Traditional DE", s=34)
    ax.scatter(plot_df["ml_rank"], y, label="Combined ML", s=34)
    ax.scatter(plot_df["rf_rank"], y, label="Random Forest", s=34)
    for _, row in plot_df.iterrows():
        yi = plot_df.index.get_loc(row.name)
        xs = [row["traditional_rank"], row["ml_rank"], row["rf_rank"]]
        ax.plot(xs, [yi, yi, yi], color="#d0d0d0", linewidth=0.8, zorder=0)

    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"])
    ax.invert_xaxis()
    ax.set_xlabel("Rank within each method; lower rank is better")
    ax.set_title("Consensus CKD-Up Biomarkers: Rank Comparison")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(figdir / "consensus_biomarker_rank_comparison.png")
    plt.close(fig)


def plot_heatmap(top: pd.DataFrame, figdir: Path) -> None:
    heat = top.set_index("gene")[
        [
            "traditional_score_scaled",
            "ml_score_scaled",
            "rf_score_scaled",
        ]
    ].copy()
    heat = heat.rename(
        columns={
            "traditional_score_scaled": "Traditional DE",
            "ml_score_scaled": "Combined ML",
            "rf_score_scaled": "Random Forest",
        }
    )
    heat = heat.iloc[::-1]

    fig, ax = plt.subplots(figsize=(6.5, max(5.5, len(heat) * 0.28)), dpi=180)
    image = ax.imshow(heat.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(np.arange(heat.shape[1]))
    ax.set_xticklabels(heat.columns, rotation=25, ha="right")
    ax.set_yticks(np.arange(heat.shape[0]))
    ax.set_yticklabels(heat.index)
    ax.set_title("Consensus CKD-Up Biomarker Evidence Heatmap")
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Scaled score")
    fig.tight_layout()
    fig.savefig(figdir / "consensus_biomarker_score_heatmap.png")
    plt.close(fig)


def plot_overlap_counts(overlap_counts: pd.DataFrame, figdir: Path) -> None:
    labels = overlap_counts["top_n"].astype(str)
    values = overlap_counts["triple_overlap_genes"].astype(int)

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=180)
    ax.bar(labels, values, color="#2f6f8f")
    ax.set_xlabel("Top N genes from each method")
    ax.set_ylabel("Triple-overlap gene count")
    ax.set_title("Traditional DE + ML + Random Forest Overlap")
    for i, value in enumerate(values):
        ax.text(i, value + max(values) * 0.02, str(value), ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(figdir / "consensus_overlap_counts.png")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    consensus, overlap_counts = load_and_merge(
        Path(args.traditional), Path(args.ml), Path(args.rf)
    )
    consensus_path = outdir / "consensus_traditional_ml_rf_biomarkers.csv"
    overlap_path = outdir / "consensus_overlap_counts.csv"
    consensus.to_csv(consensus_path, index=False)
    overlap_counts.to_csv(overlap_path, index=False)

    top = consensus.head(args.top_n).copy()
    plot_score_comparison(top, figdir)
    plot_rank_comparison(top, figdir)
    plot_heatmap(top, figdir)
    plot_overlap_counts(overlap_counts, figdir)

    print(f"Saved consensus table: {consensus_path}")
    print(f"Saved overlap counts: {overlap_path}")
    print("\nTop consensus CKD-up biomarkers:")
    print(
        top[
            [
                "gene",
                "consensus_score",
                "rank_score",
                "combined_ml_biomarker_score",
                "random_forest_importance",
                "traditional_rank",
                "ml_rank",
                "rf_rank",
                "adj_p_value",
                "diff_mean_log_expr",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
