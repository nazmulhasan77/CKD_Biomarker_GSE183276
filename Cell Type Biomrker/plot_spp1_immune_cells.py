"""
Plot SPP1+ immune cells.

SPP1+ definition:
  immune cell with SPP1 log-normalized expression > 0 in adata.raw.X.

Outputs:
  Cell Type Biomrker/results/spp1_immune_*.csv
  Cell Type Biomrker/figures/spp1_immune_*.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
from scipy import sparse

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


DEFAULT_CANDIDATES = [
    "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad",
    "../Analysis/GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad",
    "C:/Users/miaso/OneDrive/Desktop/Thesis/Analysis/GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot SPP1+ immune cells.")
    parser.add_argument("--input", default=None, help="Input h5ad path")
    parser.add_argument("--outdir", default="Cell Type Biomrker/results")
    parser.add_argument("--figdir", default="Cell Type Biomrker/figures")
    parser.add_argument("--gene", default="SPP1")
    parser.add_argument("--class-col", default="class")
    parser.add_argument("--immune-label", default="immune cells")
    parser.add_argument("--condition-col", default="condition.l1")
    parser.add_argument("--patient-col", default="patient")
    parser.add_argument("--subclass-col", default="sc.subclass.l2")
    parser.add_argument("--predicted-col", default="predicted.subclass.l3")
    return parser.parse_args()


def resolve_input(user_input: str | None) -> Path:
    candidates = [user_input] if user_input else DEFAULT_CANDIDATES
    for candidate in candidates:
        if candidate is None:
            continue
        path = Path(candidate)
        if path.exists():
            return path
    raise FileNotFoundError("Could not find h5ad file. Pass --input path.")


def get_matrix_and_genes(adata: ad.AnnData) -> tuple[object, np.ndarray]:
    if adata.raw is not None:
        matrix = adata.raw.X
        if "_index" in adata.raw.var.columns:
            genes = adata.raw.var["_index"].astype(str).to_numpy()
        else:
            genes = adata.raw.var.index.astype(str).to_numpy()
        return matrix, genes
    if "features" in adata.var.columns:
        return adata.X, adata.var["features"].astype(str).to_numpy()
    return adata.X, adata.var_names.astype(str).to_numpy()


def extract_gene_vector(matrix: object, rows: np.ndarray, gene_idx: int) -> np.ndarray:
    X = matrix[rows, [gene_idx]]
    if sparse.issparse(X):
        return np.asarray(X.toarray()).ravel()
    return np.asarray(X).ravel()


def fraction_table(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    rows = []
    for group, sub in df.groupby(group_col, observed=False):
        total = sub.shape[0]
        positive = int(sub["SPP1_positive"].sum())
        rows.append(
            {
                group_col: group,
                "n_immune_cells": total,
                "n_SPP1_positive": positive,
                "fraction_SPP1_positive": positive / total if total else 0,
                "mean_SPP1_expression": sub["SPP1_expression"].mean(),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["fraction_SPP1_positive", "n_SPP1_positive"], ascending=[False, False]
    )


def patient_table(df: pd.DataFrame, patient_col: str, condition_col: str) -> pd.DataFrame:
    rows = []
    for (patient, condition), sub in df.groupby([patient_col, condition_col], observed=False):
        total = sub.shape[0]
        positive = int(sub["SPP1_positive"].sum())
        rows.append(
            {
                patient_col: patient,
                condition_col: condition,
                "n_immune_cells": total,
                "n_SPP1_positive": positive,
                "fraction_SPP1_positive": positive / total if total else 0,
                "mean_SPP1_expression": sub["SPP1_expression"].mean(),
            }
        )
    return pd.DataFrame(rows)


def plot_umap_expression(df: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5.6), dpi=180)
    positive = df["SPP1_positive"].to_numpy()
    ax.scatter(
        df.loc[~positive, "UMAP1"],
        df.loc[~positive, "UMAP2"],
        s=2,
        c="#d0d0d0",
        alpha=0.35,
        linewidths=0,
        label="SPP1- immune",
    )
    sc = ax.scatter(
        df.loc[positive, "UMAP1"],
        df.loc[positive, "UMAP2"],
        s=6,
        c=df.loc[positive, "SPP1_expression"],
        cmap="magma",
        alpha=0.9,
        linewidths=0,
        label="SPP1+ immune",
    )
    ax.set_xlabel("ref UMAP 1")
    ax.set_ylabel("ref UMAP 2")
    ax.set_title("SPP1+ immune cells on reference UMAP")
    ax.legend(frameon=False, loc="best", markerscale=3)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("SPP1 log-normalized expression")
    fig.tight_layout()
    fig.savefig(figdir / "spp1_immune_umap_expression.png")
    plt.close(fig)


def plot_umap_condition(df: pd.DataFrame, condition_col: str, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5.6), dpi=180)
    ax.scatter(df["UMAP1"], df["UMAP2"], s=2, c="#d9d9d9", alpha=0.25, linewidths=0)
    cmap = {"CKD": "#b24a3b", "AKI": "#8060a8", "Ref": "#2f6f8f"}
    pos = df[df["SPP1_positive"]].copy()
    for cond, sub in pos.groupby(condition_col, observed=False):
        ax.scatter(
            sub["UMAP1"],
            sub["UMAP2"],
            s=7,
            c=cmap.get(str(cond), "#333333"),
            alpha=0.85,
            linewidths=0,
            label=f"{cond} SPP1+",
        )
    ax.set_xlabel("ref UMAP 1")
    ax.set_ylabel("ref UMAP 2")
    ax.set_title("SPP1+ immune cells by condition")
    ax.legend(frameon=False, loc="best", markerscale=3)
    fig.tight_layout()
    fig.savefig(figdir / "spp1_immune_umap_positive_by_condition.png")
    plt.close(fig)


def plot_fraction_by_condition(summary: pd.DataFrame, condition_col: str, figdir: Path) -> None:
    order = [x for x in ["Ref", "AKI", "CKD"] if x in set(summary[condition_col].astype(str))]
    if not order:
        order = summary[condition_col].astype(str).tolist()
    sub = summary.set_index(condition_col).loc[order].reset_index()
    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=180)
    ax.bar(sub[condition_col], sub["fraction_SPP1_positive"], color=["#2f6f8f", "#8060a8", "#b24a3b"][: len(sub)])
    ax.set_ylabel("Fraction of immune cells that are SPP1+")
    ax.set_ylim(0, max(sub["fraction_SPP1_positive"].max() * 1.25, 0.05))
    ax.set_title("SPP1+ immune cell fraction by condition")
    for i, row in sub.iterrows():
        ax.text(
            i,
            row["fraction_SPP1_positive"] + max(sub["fraction_SPP1_positive"].max() * 0.03, 0.002),
            f"{row['fraction_SPP1_positive']:.1%}\n{int(row['n_SPP1_positive'])}/{int(row['n_immune_cells'])}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    fig.tight_layout()
    fig.savefig(figdir / "spp1_immune_fraction_by_condition.png")
    plt.close(fig)


def plot_patient_fraction(patient: pd.DataFrame, condition_col: str, figdir: Path) -> None:
    order = [x for x in ["Ref", "AKI", "CKD"] if x in set(patient[condition_col].astype(str))]
    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=180)
    rng = np.random.default_rng(7)
    colors = {"Ref": "#2f6f8f", "AKI": "#8060a8", "CKD": "#b24a3b"}
    for i, cond in enumerate(order):
        vals = patient.loc[patient[condition_col].astype(str) == cond, "fraction_SPP1_positive"].to_numpy()
        jitter = rng.normal(0, 0.045, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, color=colors.get(cond, "#333333"), s=28, alpha=0.85)
        if len(vals):
            ax.hlines(np.median(vals), i - 0.22, i + 0.22, color="black", linewidth=1.4)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order)
    ax.set_ylabel("Patient-level SPP1+ immune fraction")
    ax.set_title("SPP1+ immune fraction per patient")
    fig.tight_layout()
    fig.savefig(figdir / "spp1_immune_patient_fraction_by_condition.png")
    plt.close(fig)


def plot_subclass_fraction(summary: pd.DataFrame, col: str, figdir: Path, filename: str, title: str, top_n: int = 15) -> None:
    top = summary.head(top_n).iloc[::-1]
    if top.empty:
        return
    fig, ax = plt.subplots(figsize=(7.5, max(4.5, len(top) * 0.32)), dpi=180)
    ax.barh(top[col].astype(str), top["fraction_SPP1_positive"], color="#b24a3b")
    ax.set_xlabel("Fraction SPP1+ within group")
    ax.set_title(title)
    for y, (_, row) in enumerate(top.iterrows()):
        ax.text(
            row["fraction_SPP1_positive"] + max(top["fraction_SPP1_positive"].max() * 0.02, 0.002),
            y,
            f"{int(row['n_SPP1_positive'])}/{int(row['n_immune_cells'])}",
            va="center",
            fontsize=8,
        )
    fig.tight_layout()
    fig.savefig(figdir / filename)
    plt.close(fig)


def plot_macrophage_counts(df: pd.DataFrame, subclass_col: str, condition_col: str, figdir: Path) -> None:
    mac = df[df[subclass_col].astype(str).str.contains("MAC|MON|MNP", case=False, regex=True, na=False)].copy()
    if mac.empty:
        return
    table = (
        mac.groupby([subclass_col, condition_col], observed=False)["SPP1_positive"]
        .agg(n_immune_cells="size", n_SPP1_positive="sum")
        .reset_index()
    )
    table["fraction_SPP1_positive"] = table["n_SPP1_positive"] / table["n_immune_cells"]
    table.to_csv(figdir.parent / "results" / "spp1_immune_macrophage_like_summary.csv", index=False)

    pivot = table.pivot(index=subclass_col, columns=condition_col, values="fraction_SPP1_positive").fillna(0)
    pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(7, max(4.5, pivot.shape[0] * 0.35)), dpi=180)
    x = np.arange(pivot.shape[0])
    width = 0.25
    conditions = [c for c in ["Ref", "AKI", "CKD"] if c in pivot.columns]
    colors = {"Ref": "#2f6f8f", "AKI": "#8060a8", "CKD": "#b24a3b"}
    for i, cond in enumerate(conditions):
        ax.bar(x + (i - (len(conditions) - 1) / 2) * width, pivot[cond], width=width, label=cond, color=colors.get(cond, None))
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index.astype(str), rotation=25, ha="right")
    ax.set_ylabel("Fraction SPP1+")
    ax.set_title("SPP1+ fraction in macrophage/monocyte-like immune groups")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figdir / "spp1_immune_macrophage_like_fraction.png")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    h5ad_path = resolve_input(args.input)
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path, backed="r")
    matrix, genes = get_matrix_and_genes(adata)
    gene_matches = np.flatnonzero(genes == args.gene)
    if gene_matches.size == 0:
        raise KeyError(f"{args.gene} not found in gene list")
    gene_idx = int(gene_matches[0])

    obs = adata.obs.copy()
    immune_mask = obs[args.class_col].astype(str).eq(args.immune_label).to_numpy()
    immune_rows = np.flatnonzero(immune_mask)
    expr = extract_gene_vector(matrix, immune_rows, gene_idx)
    positive = expr > 0

    coords = np.asarray(adata.obsm["X_ref.umap"])[immune_rows]
    keep_cols = [
        args.condition_col,
        args.patient_col,
        args.class_col,
        args.subclass_col,
        args.predicted_col,
        "state",
        "state.l2",
    ]
    keep_cols = [c for c in keep_cols if c in obs.columns]
    df = obs.iloc[immune_rows][keep_cols].copy()
    for col in keep_cols:
        df[col] = df[col].astype(str)
    df["UMAP1"] = coords[:, 0]
    df["UMAP2"] = coords[:, 1]
    df["SPP1_expression"] = expr
    df["SPP1_positive"] = positive
    df.to_csv(outdir / "spp1_immune_cells.csv", index=True)

    condition_summary = fraction_table(df, args.condition_col)
    patient_summary = patient_table(df, args.patient_col, args.condition_col)
    subclass_summary = fraction_table(df, args.subclass_col) if args.subclass_col in df.columns else pd.DataFrame()
    predicted_summary = fraction_table(df, args.predicted_col) if args.predicted_col in df.columns else pd.DataFrame()

    condition_summary.to_csv(outdir / "spp1_immune_condition_summary.csv", index=False)
    patient_summary.to_csv(outdir / "spp1_immune_patient_summary.csv", index=False)
    if not subclass_summary.empty:
        subclass_summary.to_csv(outdir / "spp1_immune_sc_subclass_l2_summary.csv", index=False)
    if not predicted_summary.empty:
        predicted_summary.to_csv(outdir / "spp1_immune_predicted_subclass_l3_summary.csv", index=False)

    plot_umap_expression(df, figdir)
    plot_umap_condition(df, args.condition_col, figdir)
    plot_fraction_by_condition(condition_summary, args.condition_col, figdir)
    plot_patient_fraction(patient_summary, args.condition_col, figdir)
    if not subclass_summary.empty:
        plot_subclass_fraction(
            subclass_summary,
            args.subclass_col,
            figdir,
            "spp1_immune_sc_subclass_l2_fraction.png",
            "SPP1+ immune fraction by sc.subclass.l2",
        )
        plot_macrophage_counts(df, args.subclass_col, args.condition_col, figdir)
    if not predicted_summary.empty:
        plot_subclass_fraction(
            predicted_summary,
            args.predicted_col,
            figdir,
            "spp1_immune_predicted_subclass_l3_fraction.png",
            "SPP1+ immune fraction by predicted.subclass.l3",
        )

    adata.file.close()

    print("SPP1 immune condition summary:")
    print(condition_summary.to_string(index=False))
    if not subclass_summary.empty:
        print("\nTop sc.subclass.l2 SPP1+ fractions:")
        print(subclass_summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
