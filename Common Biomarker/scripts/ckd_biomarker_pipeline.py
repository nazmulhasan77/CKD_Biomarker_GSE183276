"""
Candidate CKD biomarker discovery pipeline for GSE183276 scCv3 kidney atlas.

The input h5ad was converted from a Seurat object. In this file, adata.X is
scaled expression and can contain negative values, while adata.raw.X contains
sparse log-normalized expression. This pipeline uses adata.raw.X by default.

Main outputs:
  results/ckd_vs_ref_all_cells_de.csv
  results/ckd_vs_aki_all_cells_de.csv
  results/ckd_candidate_biomarkers_ranked.csv
  results/ckd_upregulated_candidate_biomarkers.csv
  results/ckd_downregulated_reference_loss_genes.csv
  results/subclass_l1_ckd_vs_ref_top_markers.csv
  results/patient_cell_composition.csv
  figures/*.png
"""

from __future__ import annotations

import argparse
import math
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegressionCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


DEFAULT_H5AD = "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Identify candidate CKD biomarkers from GSE183276 AnnData."
    )
    parser.add_argument("--input", default=DEFAULT_H5AD, help="Input .h5ad file")
    parser.add_argument("--outdir", default="results", help="Directory for CSV outputs")
    parser.add_argument("--figdir", default="figures", help="Directory for figures")
    parser.add_argument("--condition-col", default="condition.l1")
    parser.add_argument("--case", default="CKD")
    parser.add_argument("--control", default="Ref")
    parser.add_argument("--specificity-control", default="AKI")
    parser.add_argument("--patient-col", default="patient")
    parser.add_argument("--class-col", default="class")
    parser.add_argument("--subclass-col", default="subclass.l1")
    parser.add_argument("--chunk-size", type=int, default=1500)
    parser.add_argument("--min-pct", type=float, default=0.05)
    parser.add_argument("--min-mean", type=float, default=0.01)
    parser.add_argument("--min-cells-per-subclass-patient", type=int, default=20)
    parser.add_argument("--min-patients-per-condition", type=int, default=4)
    parser.add_argument("--top-n-subclass", type=int, default=40)
    parser.add_argument("--ml-top-n", type=int, default=500)
    parser.add_argument(
        "--skip-subclass",
        action="store_true",
        help="Skip cell-subclass-specific pseudobulk analysis.",
    )
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help="Skip UMAP, volcano, and composition figures.",
    )
    return parser.parse_args()


def safe_category(series: pd.Series) -> pd.Series:
    return series.astype("object").where(series.notna(), "NA").astype(str)


def get_expression_matrix_and_genes(adata: ad.AnnData) -> tuple[object, np.ndarray, str]:
    """Return the preferred expression matrix and gene symbols."""
    if adata.raw is not None:
        matrix = adata.raw.X
        raw_var = adata.raw.var
        if "_index" in raw_var.columns:
            genes = raw_var["_index"].astype(str).to_numpy()
        else:
            genes = raw_var.index.astype(str).to_numpy()
        return matrix, make_unique_gene_names(genes), "adata.raw.X"

    matrix = adata.X
    if "features" in adata.var.columns:
        genes = adata.var["features"].astype(str).to_numpy()
    else:
        genes = adata.var_names.astype(str).to_numpy()
    return matrix, make_unique_gene_names(genes), "adata.X"


def make_unique_gene_names(genes: np.ndarray) -> np.ndarray:
    """Make duplicate gene names unique while preserving the visible symbol."""
    seen: dict[str, int] = {}
    out: list[str] = []
    for gene in genes:
        count = seen.get(gene, 0)
        out.append(gene if count == 0 else f"{gene}__dup{count}")
        seen[gene] = count + 1
    return np.array(out, dtype=object)


def aggregate_by_groups(
    matrix: object,
    obs: pd.DataFrame,
    group_cols: list[str],
    n_genes: int,
    row_mask: np.ndarray | None = None,
    chunk_size: int = 1500,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Aggregate sparse/dense expression by metadata groups.

    Returns:
      sums: group x gene sum of log-normalized expression
      detected: group x gene count of cells with expression > 0
      group_meta: group metadata with n_cells
    """
    if row_mask is None:
        row_mask = np.ones(obs.shape[0], dtype=bool)
    selected = np.flatnonzero(row_mask)
    if selected.size == 0:
        raise ValueError("No rows selected for aggregation.")

    group_frame = obs.iloc[selected][group_cols].copy()
    for col in group_cols:
        group_frame[col] = safe_category(group_frame[col])

    key = group_frame[group_cols].agg("||".join, axis=1).to_numpy()
    unique_keys, inverse = np.unique(key, return_inverse=True)
    n_groups = unique_keys.size

    sums = np.zeros((n_groups, n_genes), dtype=np.float64)
    detected = np.zeros((n_groups, n_genes), dtype=np.uint32)
    n_cells = np.zeros(n_groups, dtype=np.int64)

    for start in range(0, selected.size, chunk_size):
        stop = min(start + chunk_size, selected.size)
        rows = selected[start:stop]
        codes = inverse[start:stop]
        X = matrix[rows, :]
        if sparse.issparse(X):
            X = X.tocsr()
        else:
            X = np.asarray(X)

        for code in np.unique(codes):
            local = np.flatnonzero(codes == code)
            Xg = X[local, :]
            n_cells[code] += local.size
            sums[code, :] += np.asarray(Xg.sum(axis=0)).ravel()
            if sparse.issparse(Xg):
                detected[code, :] += Xg.getnnz(axis=0).astype(np.uint32)
            else:
                detected[code, :] += np.count_nonzero(Xg, axis=0).astype(np.uint32)

        if start == 0 or stop == selected.size or (stop // chunk_size) % 20 == 0:
            print(f"  aggregated {stop:,}/{selected.size:,} cells")

    group_meta = pd.DataFrame(
        [key_value.split("||") for key_value in unique_keys], columns=group_cols
    )
    group_meta["n_cells"] = n_cells
    return sums, detected, group_meta


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    cleaned = np.where(np.isfinite(p_values), p_values, 1.0)
    return multipletests(cleaned, method="fdr_bh")[1]


def compare_pseudobulk_groups(
    mean_expr: np.ndarray,
    detected: np.ndarray,
    group_meta: pd.DataFrame,
    genes: np.ndarray,
    condition_col: str,
    case: str,
    control: str,
    min_pct: float,
    min_mean: float,
    min_patients_per_condition: int,
) -> pd.DataFrame:
    case_idx = group_meta[condition_col].to_numpy() == case
    control_idx = group_meta[condition_col].to_numpy() == control
    n_case = int(case_idx.sum())
    n_control = int(control_idx.sum())
    if n_case < min_patients_per_condition or n_control < min_patients_per_condition:
        raise ValueError(
            f"Need at least {min_patients_per_condition} groups per condition for "
            f"{case} vs {control}; got {n_case} and {n_control}."
        )

    case_matrix = mean_expr[case_idx, :]
    control_matrix = mean_expr[control_idx, :]
    case_mean = case_matrix.mean(axis=0)
    control_mean = control_matrix.mean(axis=0)
    diff = case_mean - control_mean

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        t_stat, p_value = stats.ttest_ind(
            case_matrix, control_matrix, axis=0, equal_var=False, nan_policy="omit"
        )

    case_sd = case_matrix.std(axis=0, ddof=1)
    control_sd = control_matrix.std(axis=0, ddof=1)
    pooled_sd = np.sqrt(
        ((n_case - 1) * case_sd**2 + (n_control - 1) * control_sd**2)
        / max(n_case + n_control - 2, 1)
    )
    cohen_d = np.divide(
        diff,
        pooled_sd,
        out=np.zeros_like(diff, dtype=np.float64),
        where=pooled_sd > 0,
    )

    case_cells = group_meta.loc[case_idx, "n_cells"].astype(int).to_numpy()
    control_cells = group_meta.loc[control_idx, "n_cells"].astype(int).to_numpy()
    pct_case = detected[case_idx, :].sum(axis=0) / max(case_cells.sum(), 1)
    pct_control = detected[control_idx, :].sum(axis=0) / max(control_cells.sum(), 1)

    adj_p_value = benjamini_hochberg(p_value)
    max_pct = np.maximum(pct_case, pct_control)
    max_mean = np.maximum(case_mean, control_mean)
    passes_expression = (max_pct >= min_pct) & (max_mean >= min_mean)
    rank_score = (
        -np.log10(np.clip(adj_p_value, 1e-300, 1.0))
        * np.abs(cohen_d)
        * np.abs(diff)
    )
    rank_score = np.where(passes_expression, rank_score, 0.0)

    out = pd.DataFrame(
        {
            "gene": genes,
            "mean_expr_case": case_mean,
            "mean_expr_control": control_mean,
            "diff_mean_log_expr": diff,
            "cohen_d": cohen_d,
            "t_stat": t_stat,
            "p_value": p_value,
            "adj_p_value": adj_p_value,
            "pct_expr_case": pct_case,
            "pct_expr_control": pct_control,
            "pct_expr_diff": pct_case - pct_control,
            "n_case_groups": n_case,
            "n_control_groups": n_control,
            "passes_expression_filter": passes_expression,
            "direction": np.where(diff >= 0, "up_in_case", "down_in_case"),
            "rank_score": rank_score,
        }
    )
    out = out.sort_values(
        ["passes_expression_filter", "rank_score", "adj_p_value"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    return out


def run_elasticnet_marker_model(
    mean_expr: np.ndarray,
    group_meta: pd.DataFrame,
    de_table: pd.DataFrame,
    genes: np.ndarray,
    condition_col: str,
    case: str,
    control: str,
    top_n: int,
) -> pd.DataFrame:
    """Patient-level logistic model for an additional feature-importance view."""
    rows = group_meta[condition_col].isin([case, control]).to_numpy()
    y = (group_meta.loc[rows, condition_col].to_numpy() == case).astype(int)
    selected_genes = (
        de_table.loc[de_table["passes_expression_filter"], "gene"].head(top_n).tolist()
    )
    if len(selected_genes) < 5:
        return pd.DataFrame()

    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    gene_idx = np.array([gene_to_idx[gene] for gene in selected_genes], dtype=int)
    X = mean_expr[rows, :][:, gene_idx]

    class_counts = np.bincount(y)
    cv_folds = int(min(5, class_counts.min()))
    if cv_folds < 3:
        return pd.DataFrame()

    model = make_pipeline(
        StandardScaler(),
        LogisticRegressionCV(
            Cs=10,
            cv=cv_folds,
            penalty="elasticnet",
            solver="saga",
            l1_ratios=[0.25, 0.5, 0.75, 1.0],
            scoring="balanced_accuracy",
            class_weight="balanced",
            max_iter=20000,
            n_jobs=1,
            random_state=7,
        ),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        model.fit(X, y)

    clf = model.named_steps["logisticregressioncv"]
    coef = clf.coef_.ravel()
    importance = pd.DataFrame(
        {
            "gene": selected_genes,
            "elasticnet_coefficient": coef,
            "abs_elasticnet_coefficient": np.abs(coef),
            "selected_by_model": np.abs(coef) > 1e-9,
            "best_C": float(clf.C_[0]),
            "best_l1_ratio": float(clf.l1_ratio_[0]),
        }
    )
    importance = importance.merge(
        de_table[
            [
                "gene",
                "rank_score",
                "adj_p_value",
                "diff_mean_log_expr",
                "cohen_d",
                "pct_expr_case",
                "pct_expr_control",
            ]
        ],
        on="gene",
        how="left",
    )
    return importance.sort_values(
        ["selected_by_model", "abs_elasticnet_coefficient", "rank_score"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


def save_metadata_summaries(
    obs: pd.DataFrame,
    outdir: Path,
    condition_col: str,
    patient_col: str,
    class_col: str,
    subclass_col: str,
) -> pd.DataFrame:
    patient_summary = (
        obs.groupby([patient_col, condition_col], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values([condition_col, patient_col])
    )
    patient_summary.to_csv(outdir / "patient_summary.csv", index=False)

    class_counts = (
        obs.groupby([patient_col, condition_col, class_col], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    totals = class_counts.groupby([patient_col, condition_col], observed=False)[
        "n_cells"
    ].transform("sum")
    class_counts["fraction"] = class_counts["n_cells"] / totals
    class_counts.to_csv(outdir / "patient_cell_composition.csv", index=False)

    subclass_counts = (
        obs.groupby([patient_col, condition_col, subclass_col], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    subclass_counts.to_csv(outdir / "patient_subclass_counts.csv", index=False)
    return class_counts


def plot_umap(
    adata: ad.AnnData,
    obs: pd.DataFrame,
    figdir: Path,
    condition_col: str,
    class_col: str,
) -> None:
    if "X_ref.umap" not in adata.obsm:
        return
    umap = np.asarray(adata.obsm["X_ref.umap"])
    plot_categorical_umap(umap, obs[condition_col], figdir / "umap_condition_l1.png")
    plot_categorical_umap(umap, obs[class_col], figdir / "umap_class.png")


def plot_categorical_umap(coords: np.ndarray, labels: pd.Series, path: Path) -> None:
    cats = pd.Categorical(labels)
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(8, 6), dpi=180)
    for idx, cat in enumerate(cats.categories):
        mask = cats == cat
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=1,
            alpha=0.5,
            color=cmap(idx % 20),
            label=str(cat),
            linewidths=0,
        )
    ax.set_xlabel("ref UMAP 1")
    ax.set_ylabel("ref UMAP 2")
    ax.set_title(path.stem.replace("_", " "))
    ax.legend(markerscale=6, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_composition(
    class_counts: pd.DataFrame, figdir: Path, condition_col: str, class_col: str
) -> None:
    summary = (
        class_counts.groupby([condition_col, class_col], observed=False)["fraction"]
        .mean()
        .reset_index()
    )
    pivot = summary.pivot(
        index=condition_col, columns=class_col, values="fraction"
    ).fillna(0)
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=180)
    bottom = np.zeros(pivot.shape[0])
    colors = plt.get_cmap("Set2")(np.linspace(0, 1, pivot.shape[1]))
    for i, col in enumerate(pivot.columns):
        ax.bar(pivot.index, pivot[col], bottom=bottom, label=col, color=colors[i])
        bottom += pivot[col].to_numpy()
    ax.set_ylabel("Mean patient cell fraction")
    ax.set_xlabel("")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(figdir / "cell_class_composition_by_condition.png")
    plt.close(fig)


def plot_volcano(de_table: pd.DataFrame, figdir: Path, case: str, control: str) -> None:
    df = de_table.copy()
    df["neg_log10_fdr"] = -np.log10(np.clip(df["adj_p_value"], 1e-300, 1.0))
    sig_up = (
        (df["adj_p_value"] < 0.05)
        & (df["diff_mean_log_expr"] > 0)
        & df["passes_expression_filter"]
    )
    sig_down = (
        (df["adj_p_value"] < 0.05)
        & (df["diff_mean_log_expr"] < 0)
        & df["passes_expression_filter"]
    )
    colors = np.where(sig_up, "#c0392b", np.where(sig_down, "#2878b5", "#b8b8b8"))
    fig, ax = plt.subplots(figsize=(7, 5), dpi=180)
    ax.scatter(
        df["diff_mean_log_expr"],
        df["neg_log10_fdr"],
        s=5,
        c=colors,
        alpha=0.65,
        linewidths=0,
    )
    ax.axhline(-math.log10(0.05), color="#555555", linewidth=0.8, linestyle=":")
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_xlabel(f"Mean log-expression difference ({case} - {control})")
    ax.set_ylabel("-log10 FDR")
    ax.set_title(f"{case} vs {control} patient-level pseudobulk")

    label_df = pd.concat(
        [
            df.loc[sig_up].head(10),
            df.loc[sig_down].head(10),
        ],
        axis=0,
    )
    for _, row in label_df.iterrows():
        ax.text(
            row["diff_mean_log_expr"],
            row["neg_log10_fdr"],
            row["gene"],
            fontsize=6,
        )
    fig.tight_layout()
    fig.savefig(figdir / "volcano_ckd_vs_ref_all_cells.png")
    plt.close(fig)


def plot_top_biomarkers(candidate_table: pd.DataFrame, figdir: Path) -> None:
    for direction, title, path, color in [
        (
            "up",
            "Top CKD-upregulated candidate biomarkers",
            "top20_ckd_upregulated_candidate_biomarkers.png",
            "#c0392b",
        ),
        (
            "down",
            "Top CKD-downregulated/reference-loss genes",
            "top20_ckd_downregulated_reference_loss_genes.png",
            "#2878b5",
        ),
    ]:
        if direction == "up":
            top = candidate_table[candidate_table["diff_mean_log_expr"] > 0].head(20)
        else:
            top = candidate_table[candidate_table["diff_mean_log_expr"] < 0].head(20)
        if top.empty:
            continue
        top = top.iloc[::-1]
        fig, ax = plt.subplots(figsize=(7, 6), dpi=180)
        ax.barh(top["gene"], top["rank_score"], color=color)
        ax.set_xlabel("Composite biomarker rank score")
        ax.set_ylabel("")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(figdir / path)
        plt.close(fig)


def add_specificity_columns(
    ckd_ref: pd.DataFrame,
    ckd_aki: pd.DataFrame,
    case: str,
    specificity_control: str,
) -> pd.DataFrame:
    specificity = ckd_aki[
        [
            "gene",
            "diff_mean_log_expr",
            "adj_p_value",
            "pct_expr_case",
            "pct_expr_control",
        ]
    ].rename(
        columns={
            "diff_mean_log_expr": f"diff_mean_log_expr_{case}_vs_{specificity_control}",
            "adj_p_value": f"adj_p_value_{case}_vs_{specificity_control}",
            "pct_expr_case": f"pct_expr_{case}",
            "pct_expr_control": f"pct_expr_{specificity_control}",
        }
    )
    out = ckd_ref.merge(specificity, on="gene", how="left")
    out["ckd_specificity_score"] = np.where(
        out[f"diff_mean_log_expr_{case}_vs_{specificity_control}"].fillna(0) > 0,
        -np.log10(
            np.clip(
                out[f"adj_p_value_{case}_vs_{specificity_control}"].fillna(1.0),
                1e-300,
                1.0,
            )
        ),
        0.0,
    )
    out["combined_biomarker_score"] = out["rank_score"] * (
        1.0 + out["ckd_specificity_score"]
    )
    return out.sort_values(
        [
            "passes_expression_filter",
            "combined_biomarker_score",
            "rank_score",
            "adj_p_value",
        ],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)


def save_directional_candidate_tables(
    candidates: pd.DataFrame,
    outdir: Path,
    case: str,
    specificity_control: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    up = candidates[
        (candidates["diff_mean_log_expr"] > 0)
        & candidates["passes_expression_filter"]
    ].copy()
    down = candidates[
        (candidates["diff_mean_log_expr"] < 0)
        & candidates["passes_expression_filter"]
    ].copy()

    specificity_diff_col = f"diff_mean_log_expr_{case}_vs_{specificity_control}"
    specificity_fdr_col = f"adj_p_value_{case}_vs_{specificity_control}"
    strict_specific = up[
        (up[specificity_diff_col] > 0) & (up[specificity_fdr_col] < 0.05)
    ].copy()

    up.to_csv(outdir / "ckd_upregulated_candidate_biomarkers.csv", index=False)
    down.to_csv(outdir / "ckd_downregulated_reference_loss_genes.csv", index=False)
    strict_specific.to_csv(
        outdir / "ckd_upregulated_and_aki_specific_biomarkers.csv", index=False
    )
    return up, down, strict_specific


def subclass_marker_table(
    mean_expr: np.ndarray,
    detected: np.ndarray,
    meta: pd.DataFrame,
    genes: np.ndarray,
    condition_col: str,
    subclass_col: str,
    case: str,
    control: str,
    min_pct: float,
    min_mean: float,
    min_cells_per_subclass_patient: int,
    min_patients_per_condition: int,
    top_n: int,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    subclasses = sorted(meta[subclass_col].unique())
    for subclass in subclasses:
        sub_mask = (meta[subclass_col] == subclass) & (
            meta["n_cells"] >= min_cells_per_subclass_patient
        )
        sub_meta = meta.loc[sub_mask].reset_index(drop=True)
        if sub_meta.empty:
            continue
        sub_mean = mean_expr[sub_mask.to_numpy(), :]
        sub_detected = detected[sub_mask.to_numpy(), :]
        counts = sub_meta[condition_col].value_counts()
        if counts.get(case, 0) < min_patients_per_condition:
            continue
        if counts.get(control, 0) < min_patients_per_condition:
            continue
        try:
            de = compare_pseudobulk_groups(
                sub_mean,
                sub_detected,
                sub_meta,
                genes,
                condition_col,
                case,
                control,
                min_pct,
                min_mean,
                min_patients_per_condition,
            )
        except ValueError:
            continue
        de.insert(0, subclass_col, subclass)

        top_overall = de.head(top_n).copy()
        top_overall.insert(1, "marker_set", "top_overall")

        top_up = de[de["diff_mean_log_expr"] > 0].head(top_n).copy()
        top_up.insert(1, "marker_set", "up_in_CKD")

        top_down = de[de["diff_mean_log_expr"] < 0].head(top_n).copy()
        top_down.insert(1, "marker_set", "down_in_CKD")

        rows.append(pd.concat([top_overall, top_up, top_down], axis=0))
        print(
            f"  subclass {subclass}: kept top {min(top_n, de.shape[0])} overall, "
            f"{top_up.shape[0]} up, {top_down.shape[0]} down markers"
        )
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, axis=0, ignore_index=True)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.input}")
    adata = ad.read_h5ad(args.input, backed="r")
    obs = adata.obs.copy()
    for col in [
        args.condition_col,
        args.patient_col,
        args.class_col,
        args.subclass_col,
    ]:
        if col not in obs.columns:
            raise KeyError(f"Missing required obs column: {col}")
        obs[col] = safe_category(obs[col])

    matrix, genes, matrix_name = get_expression_matrix_and_genes(adata)
    print(f"Using {matrix_name} with {len(genes):,} genes.")

    print("Saving metadata summaries.")
    class_counts = save_metadata_summaries(
        obs,
        outdir,
        args.condition_col,
        args.patient_col,
        args.class_col,
        args.subclass_col,
    )

    conditions_for_all = [args.case, args.control, args.specificity_control]
    all_mask = obs[args.condition_col].isin(conditions_for_all).to_numpy()
    print("Aggregating all cells by patient.")
    all_sums, all_detected, all_meta = aggregate_by_groups(
        matrix,
        obs,
        [args.patient_col, args.condition_col],
        len(genes),
        row_mask=all_mask,
        chunk_size=args.chunk_size,
    )
    all_mean = all_sums / np.maximum(all_meta["n_cells"].to_numpy()[:, None], 1)
    all_meta.to_csv(outdir / "pseudobulk_all_cells_metadata.csv", index=False)

    print(f"Testing {args.case} vs {args.control}.")
    ckd_ref = compare_pseudobulk_groups(
        all_mean,
        all_detected,
        all_meta,
        genes,
        args.condition_col,
        args.case,
        args.control,
        args.min_pct,
        args.min_mean,
        args.min_patients_per_condition,
    )
    ckd_ref.to_csv(outdir / "ckd_vs_ref_all_cells_de.csv", index=False)

    print(f"Testing {args.case} vs {args.specificity_control}.")
    ckd_aki = compare_pseudobulk_groups(
        all_mean,
        all_detected,
        all_meta,
        genes,
        args.condition_col,
        args.case,
        args.specificity_control,
        args.min_pct,
        args.min_mean,
        args.min_patients_per_condition,
    )
    ckd_aki.to_csv(outdir / "ckd_vs_aki_all_cells_de.csv", index=False)

    candidates = add_specificity_columns(
        ckd_ref, ckd_aki, args.case, args.specificity_control
    )
    candidates.to_csv(outdir / "ckd_candidate_biomarkers_ranked.csv", index=False)
    up_candidates, down_candidates, strict_specific = save_directional_candidate_tables(
        candidates, outdir, args.case, args.specificity_control
    )

    print("Running patient-level elastic-net feature importance.")
    ml_importance = run_elasticnet_marker_model(
        all_mean,
        all_meta,
        ckd_ref,
        genes,
        args.condition_col,
        args.case,
        args.control,
        args.ml_top_n,
    )
    if not ml_importance.empty:
        ml_importance.to_csv(outdir / "elasticnet_gene_importance.csv", index=False)

    if not args.skip_subclass:
        print("Aggregating by patient and subclass.")
        subclass_mask = obs[args.condition_col].isin([args.case, args.control]).to_numpy()
        sub_sums, sub_detected, sub_meta = aggregate_by_groups(
            matrix,
            obs,
            [args.subclass_col, args.patient_col, args.condition_col],
            len(genes),
            row_mask=subclass_mask,
            chunk_size=args.chunk_size,
        )
        sub_mean = sub_sums / np.maximum(sub_meta["n_cells"].to_numpy()[:, None], 1)
        sub_meta.to_csv(outdir / "pseudobulk_subclass_l1_metadata.csv", index=False)
        sub_markers = subclass_marker_table(
            sub_mean,
            sub_detected,
            sub_meta,
            genes,
            args.condition_col,
            args.subclass_col,
            args.case,
            args.control,
            args.min_pct,
            args.min_mean,
            args.min_cells_per_subclass_patient,
            args.min_patients_per_condition,
            args.top_n_subclass,
        )
        if not sub_markers.empty:
            sub_markers.to_csv(
                outdir / "subclass_l1_ckd_vs_ref_top_markers.csv", index=False
            )

    if not args.skip_figures:
        print("Creating figures.")
        plot_umap(adata, obs, figdir, args.condition_col, args.class_col)
        plot_composition(class_counts, figdir, args.condition_col, args.class_col)
        plot_volcano(ckd_ref, figdir, args.case, args.control)
        plot_top_biomarkers(candidates, figdir)

    adata.file.close()
    print("Done. Top CKD candidate biomarkers:")
    if strict_specific.empty:
        print(
            "No all-cell CKD-up gene passed strict CKD-vs-Ref FDR < 0.05 and "
            f"{args.case}-vs-{args.specificity_control} FDR < 0.05. "
            "Use the CKD-vs-Ref upregulated table as candidate biomarkers and "
            "validate CKD specificity externally."
        )
    print("\nTop CKD-upregulated candidates:")
    print(
        up_candidates[
            [
                "gene",
                "combined_biomarker_score",
                "rank_score",
                "adj_p_value",
                "diff_mean_log_expr",
                "cohen_d",
                "pct_expr_case",
                "pct_expr_control",
            ]
        ]
        .head(25)
        .to_string(index=False)
    )
    print("\nTop CKD-downregulated/reference-loss genes:")
    print(
        down_candidates[
            [
                "gene",
                "combined_biomarker_score",
                "rank_score",
                "adj_p_value",
                "diff_mean_log_expr",
                "cohen_d",
                "pct_expr_case",
                "pct_expr_control",
            ]
        ]
        .head(25)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
