"""
Cell-type-specific CKD biomarker discovery pipeline.

For each major cell class, this script performs:
  1. Patient-level pseudobulk aggregation within that cell class
  2. Traditional CKD vs Ref differential expression
  3. Random Forest feature importance
  4. KNN single-gene cross-validation importance
  5. Consensus biomarker ranking and plots

Outputs are written under:
  Cell Type Biomrker/results
  Cell Type Biomrker/figures
  Cell Type Biomrker/data
"""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


DEFAULT_H5AD = "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad"
DEFAULT_CLASSES = [
    "epithelial cells",
    "immune cells",
    "endothelial cells",
    "stroma cells",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cell-type-specific CKD biomarker analysis."
    )
    parser.add_argument("--input", default=DEFAULT_H5AD)
    parser.add_argument("--outdir", default="Cell Type Biomrker/results")
    parser.add_argument("--figdir", default="Cell Type Biomrker/figures")
    parser.add_argument("--datadir", default="Cell Type Biomrker/data")
    parser.add_argument("--condition-col", default="condition.l1")
    parser.add_argument("--patient-col", default="patient")
    parser.add_argument("--class-col", default="class")
    parser.add_argument("--case", default="CKD")
    parser.add_argument("--control", default="Ref")
    parser.add_argument("--chunk-size", type=int, default=2000)
    parser.add_argument("--min-pct", type=float, default=0.03)
    parser.add_argument("--min-mean", type=float, default=0.005)
    parser.add_argument("--min-cells-per-patient", type=int, default=10)
    parser.add_argument("--top-n-filter", type=int, default=300)
    parser.add_argument("--rf-trees", type=int, default=1000)
    parser.add_argument("--permutation-repeats", type=int, default=30)
    parser.add_argument("--random-state", type=int, default=7)
    parser.add_argument("--plot-top-n", type=int, default=20)
    return parser.parse_args()


def slugify(text: str) -> str:
    text = text.lower().strip()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def make_unique_gene_names(genes: np.ndarray) -> np.ndarray:
    seen: dict[str, int] = {}
    out: list[str] = []
    for gene in genes.astype(str):
        count = seen.get(gene, 0)
        out.append(gene if count == 0 else f"{gene}__dup{count}")
        seen[gene] = count + 1
    return np.array(out, dtype=object)


def get_expression_matrix_and_genes(adata: ad.AnnData) -> tuple[object, np.ndarray]:
    if adata.raw is not None:
        matrix = adata.raw.X
        raw_var = adata.raw.var
        if "_index" in raw_var.columns:
            genes = raw_var["_index"].astype(str).to_numpy()
        else:
            genes = raw_var.index.astype(str).to_numpy()
        return matrix, make_unique_gene_names(genes)

    matrix = adata.X
    if "features" in adata.var.columns:
        genes = adata.var["features"].astype(str).to_numpy()
    else:
        genes = adata.var_names.astype(str).to_numpy()
    return matrix, make_unique_gene_names(genes)


def aggregate_class_patient(
    matrix: object,
    obs: pd.DataFrame,
    genes: np.ndarray,
    cell_class: str,
    patient_col: str,
    condition_col: str,
    class_col: str,
    case: str,
    control: str,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    row_mask = (
        (obs[class_col].astype(str).to_numpy() == cell_class)
        & obs[condition_col].isin([case, control]).to_numpy()
    )
    selected = np.flatnonzero(row_mask)
    if selected.size == 0:
        raise ValueError(f"No rows selected for class: {cell_class}")

    group_obs = obs.iloc[selected][[patient_col, condition_col, class_col]].copy()
    for col in [patient_col, condition_col, class_col]:
        group_obs[col] = group_obs[col].astype(str)
    keys = (group_obs[patient_col] + "||" + group_obs[condition_col]).to_numpy()
    unique_keys, inverse = np.unique(keys, return_inverse=True)

    n_groups = unique_keys.size
    n_genes = len(genes)
    sums = np.zeros((n_groups, n_genes), dtype=np.float64)
    detected = np.zeros((n_groups, n_genes), dtype=np.uint32)
    n_cells = np.zeros(n_groups, dtype=np.int64)

    print(f"\n{cell_class}: aggregating {selected.size:,} cells into {n_groups} patients")
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
            sums[code, :] += np.asarray(Xg.sum(axis=0)).ravel()
            if sparse.issparse(Xg):
                detected[code, :] += Xg.getnnz(axis=0).astype(np.uint32)
            else:
                detected[code, :] += np.count_nonzero(Xg, axis=0).astype(np.uint32)
            n_cells[code] += local.size

        if start == 0 or stop == selected.size or (stop // chunk_size) % 20 == 0:
            print(f"  processed {stop:,}/{selected.size:,} cells")

    meta = pd.DataFrame(
        [key.split("||") for key in unique_keys],
        columns=[patient_col, condition_col],
    )
    meta[class_col] = cell_class
    meta["n_cells"] = n_cells
    return sums, detected, meta


def compare_de(
    mean_expr: np.ndarray,
    detected: np.ndarray,
    meta: pd.DataFrame,
    genes: np.ndarray,
    condition_col: str,
    case: str,
    control: str,
    min_pct: float,
    min_mean: float,
) -> pd.DataFrame:
    case_idx = meta[condition_col].to_numpy() == case
    control_idx = meta[condition_col].to_numpy() == control
    case_x = mean_expr[case_idx]
    control_x = mean_expr[control_idx]
    n_case = case_x.shape[0]
    n_control = control_x.shape[0]

    case_mean = case_x.mean(axis=0)
    control_mean = control_x.mean(axis=0)
    diff = case_mean - control_mean
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        t_stat, p_value = stats.ttest_ind(
            case_x, control_x, axis=0, equal_var=False, nan_policy="omit"
        )
    p_value = np.nan_to_num(p_value, nan=1.0, posinf=1.0, neginf=1.0)
    adj_p = multipletests(p_value, method="fdr_bh")[1]

    case_sd = case_x.std(axis=0, ddof=1)
    control_sd = control_x.std(axis=0, ddof=1)
    pooled_sd = np.sqrt(
        ((n_case - 1) * case_sd**2 + (n_control - 1) * control_sd**2)
        / max(n_case + n_control - 2, 1)
    )
    cohen_d = np.divide(diff, pooled_sd, out=np.zeros_like(diff), where=pooled_sd > 0)

    case_cells = meta.loc[case_idx, "n_cells"].astype(int).sum()
    control_cells = meta.loc[control_idx, "n_cells"].astype(int).sum()
    pct_case = detected[case_idx].sum(axis=0) / max(case_cells, 1)
    pct_control = detected[control_idx].sum(axis=0) / max(control_cells, 1)

    passes = (np.maximum(pct_case, pct_control) >= min_pct) & (
        np.maximum(case_mean, control_mean) >= min_mean
    )
    rank_score = (
        -np.log10(np.clip(adj_p, 1e-300, 1.0))
        * np.abs(cohen_d)
        * np.abs(diff)
    )
    rank_score = np.where(passes, rank_score, 0.0)

    out = pd.DataFrame(
        {
            "gene": genes,
            "mean_expr_case": case_mean,
            "mean_expr_control": control_mean,
            "diff_mean_log_expr": diff,
            "direction": np.where(diff >= 0, "up_in_CKD", "down_in_CKD"),
            "cohen_d": cohen_d,
            "t_stat": t_stat,
            "p_value": p_value,
            "adj_p_value": adj_p,
            "pct_expr_case": pct_case,
            "pct_expr_control": pct_control,
            "pct_expr_diff": pct_case - pct_control,
            "n_case_patients": n_case,
            "n_control_patients": n_control,
            "passes_expression_filter": passes,
            "rank_score": rank_score,
        }
    )
    return out.sort_values(
        ["passes_expression_filter", "rank_score", "adj_p_value"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def safe_f_classif(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y)
    classes = np.unique(y)
    n_samples = X.shape[0]
    n_classes = classes.size
    overall_mean = X.mean(axis=0)
    ss_between = np.zeros(X.shape[1], dtype=np.float64)
    ss_within = np.zeros(X.shape[1], dtype=np.float64)
    for cls in classes:
        Xc = X[y == cls]
        class_mean = Xc.mean(axis=0)
        ss_between += Xc.shape[0] * (class_mean - overall_mean) ** 2
        ss_within += ((Xc - class_mean) ** 2).sum(axis=0)
    df_between = max(n_classes - 1, 1)
    df_within = max(n_samples - n_classes, 1)
    scores = np.divide(
        ss_between / df_between,
        ss_within / df_within,
        out=np.zeros_like(ss_between),
        where=ss_within > 0,
    )
    p_values = stats.f.sf(scores, df_between, df_within)
    scores = np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)
    p_values = np.nan_to_num(p_values, nan=1.0, posinf=1.0, neginf=1.0)
    return scores, p_values


def make_cv(y: np.ndarray, random_state: int) -> StratifiedKFold:
    n_splits = int(min(5, np.bincount(y).min()))
    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)


def model_cv_metrics(
    X: np.ndarray,
    y: np.ndarray,
    top_n: int,
    rf_trees: int,
    random_state: int,
) -> pd.DataFrame:
    k = min(top_n, X.shape[1])
    cv = make_cv(y, random_state)
    models = {
        "RandomForest": Pipeline(
            [
                ("select", SelectKBest(score_func=safe_f_classif, k=k)),
                (
                    "model",
                    RandomForestClassifier(
                        n_estimators=rf_trees,
                        random_state=random_state,
                        class_weight="balanced",
                        max_features="sqrt",
                        n_jobs=-1,
                    ),
                ),
            ]
        ),
        "KNN": Pipeline(
            [
                ("select", SelectKBest(score_func=safe_f_classif, k=k)),
                ("scale", StandardScaler()),
                ("model", KNeighborsClassifier(n_neighbors=5, weights="distance")),
            ]
        ),
    }
    scoring = {
        "accuracy": "accuracy",
        "balanced_accuracy": "balanced_accuracy",
        "f1": "f1",
        "roc_auc": "roc_auc",
    }
    rows: list[dict[str, object]] = []
    for model_name, model in models.items():
        scores = cross_validate(
            model,
            X,
            y,
            cv=cv,
            scoring=scoring,
            n_jobs=-1,
            error_score=np.nan,
        )
        for metric in scoring:
            vals = scores[f"test_{metric}"]
            rows.append(
                {
                    "model": model_name,
                    "metric": metric,
                    "mean": float(np.nanmean(vals)),
                    "std": float(np.nanstd(vals, ddof=1)),
                    "fold_values": ";".join(f"{v:.4f}" for v in vals),
                }
            )
    return pd.DataFrame(rows)


def fit_ml(
    mean_expr: np.ndarray,
    meta: pd.DataFrame,
    genes: np.ndarray,
    condition_col: str,
    case: str,
    control: str,
    top_n_filter: int,
    rf_trees: int,
    permutation_repeats: int,
    random_state: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = meta[condition_col].isin([case, control]).to_numpy()
    X = mean_expr[rows]
    y = (meta.loc[rows, condition_col].to_numpy() == case).astype(int)
    keep = X.var(axis=0) > 0
    X = X[:, keep]
    kept_genes = genes[keep]

    cv_metrics = model_cv_metrics(X, y, top_n_filter, rf_trees, random_state)

    k = min(top_n_filter, X.shape[1])
    selector = SelectKBest(score_func=safe_f_classif, k=k)
    selector.fit(X, y)
    idx = selector.get_support(indices=True)
    X_sel = X[:, idx]
    selected_genes = kept_genes[idx]

    rf = RandomForestClassifier(
        n_estimators=rf_trees,
        random_state=random_state,
        class_weight="balanced",
        max_features="sqrt",
        n_jobs=-1,
    )
    rf.fit(X_sel, y)
    rf_importance = pd.DataFrame(
        {"gene": selected_genes, "random_forest_importance": rf.feature_importances_}
    ).sort_values("random_forest_importance", ascending=False)

    knn = Pipeline(
        [
            ("scale", StandardScaler()),
            ("model", KNeighborsClassifier(n_neighbors=5, weights="distance")),
        ]
    )
    knn.fit(X_sel, y)
    perm = permutation_importance(
        knn,
        X_sel,
        y,
        scoring="balanced_accuracy",
        n_repeats=permutation_repeats,
        random_state=random_state,
        n_jobs=-1,
    )
    knn_perm = pd.DataFrame(
        {
            "gene": selected_genes,
            "knn_permutation_importance_mean": perm.importances_mean,
            "knn_permutation_importance_std": perm.importances_std,
        }
    ).sort_values("knn_permutation_importance_mean", ascending=False)

    cv = make_cv(y, random_state)
    knn_single_rows: list[dict[str, object]] = []
    single_gene_model = Pipeline(
        [
            ("scale", StandardScaler()),
            ("model", KNeighborsClassifier(n_neighbors=5, weights="distance")),
        ]
    )
    scoring = {
        "accuracy": "accuracy",
        "balanced_accuracy": "balanced_accuracy",
        "roc_auc": "roc_auc",
    }
    for col_idx, gene in enumerate(selected_genes):
        scores = cross_validate(
            single_gene_model,
            X_sel[:, [col_idx]],
            y,
            cv=cv,
            scoring=scoring,
            n_jobs=1,
            error_score=np.nan,
        )
        balanced = float(np.nanmean(scores["test_balanced_accuracy"]))
        knn_single_rows.append(
            {
                "gene": gene,
                "knn_single_gene_accuracy_mean": float(np.nanmean(scores["test_accuracy"])),
                "knn_single_gene_balanced_accuracy_mean": balanced,
                "knn_single_gene_roc_auc_mean": float(np.nanmean(scores["test_roc_auc"])),
                "knn_single_gene_importance_score": max(balanced - 0.5, 0.0),
            }
        )
    knn_single = pd.DataFrame(knn_single_rows).sort_values(
        "knn_single_gene_importance_score", ascending=False
    )

    final_fit = pd.DataFrame(
        [
            summarize_final_model("RandomForest_final_fit", rf, X_sel, y),
            summarize_final_model("KNN_final_fit", knn, X_sel, y),
        ]
    )
    return (
        rf_importance.reset_index(drop=True),
        knn_perm.reset_index(drop=True),
        knn_single.reset_index(drop=True),
        cv_metrics,
        final_fit,
    )


def summarize_final_model(model_name: str, model: object, X: np.ndarray, y: np.ndarray) -> dict[str, object]:
    pred = model.predict(X)
    prob = model.predict_proba(X)[:, 1]
    return {
        "model": model_name,
        "accuracy": accuracy_score(y, pred),
        "balanced_accuracy": balanced_accuracy_score(y, pred),
        "f1": f1_score(y, pred),
        "roc_auc": roc_auc_score(y, prob),
    }


def minmax(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    lo = values.min()
    hi = values.max()
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        return pd.Series(np.ones(len(values)), index=values.index)
    return (values - lo) / (hi - lo)


def combine_ml(de: pd.DataFrame, rf: pd.DataFrame, knn_perm: pd.DataFrame, knn_single: pd.DataFrame) -> pd.DataFrame:
    stats_cols = [
        "gene",
        "diff_mean_log_expr",
        "direction",
        "cohen_d",
        "adj_p_value",
        "rank_score",
        "pct_expr_case",
        "pct_expr_control",
    ]
    out = de[stats_cols].merge(rf, on="gene", how="inner")
    out = out.merge(knn_perm, on="gene", how="inner")
    out = out.merge(knn_single, on="gene", how="inner")
    out["traditional_score_scaled"] = minmax(out["rank_score"])
    out["rf_score_scaled"] = minmax(out["random_forest_importance"])
    out["knn_score_scaled"] = minmax(out["knn_single_gene_importance_score"])
    out["combined_ml_biomarker_score"] = (
        0.40 * out["rf_score_scaled"]
        + 0.40 * out["knn_score_scaled"]
        + 0.20 * out["traditional_score_scaled"]
    )
    return out.sort_values("combined_ml_biomarker_score", ascending=False).reset_index(drop=True)


def consensus_table(de_up: pd.DataFrame, ml_up: pd.DataFrame, rf: pd.DataFrame) -> pd.DataFrame:
    de_ranked = de_up.sort_values(["rank_score", "adj_p_value"], ascending=[False, True]).copy()
    ml_ranked = ml_up.sort_values("combined_ml_biomarker_score", ascending=False).copy()
    rf_ranked = rf.sort_values("random_forest_importance", ascending=False).copy()
    de_ranked["traditional_rank"] = np.arange(1, len(de_ranked) + 1)
    ml_ranked["ml_rank"] = np.arange(1, len(ml_ranked) + 1)
    rf_ranked["rf_rank"] = np.arange(1, len(rf_ranked) + 1)

    con = de_ranked[
        [
            "gene",
            "rank_score",
            "adj_p_value",
            "diff_mean_log_expr",
            "cohen_d",
            "pct_expr_case",
            "pct_expr_control",
            "traditional_rank",
        ]
    ].merge(
        ml_ranked[
            [
                "gene",
                "combined_ml_biomarker_score",
                "knn_single_gene_balanced_accuracy_mean",
                "knn_single_gene_roc_auc_mean",
                "ml_rank",
            ]
        ],
        on="gene",
        how="inner",
    )
    con = con.merge(
        rf_ranked[["gene", "random_forest_importance", "rf_rank"]],
        on="gene",
        how="inner",
    )
    con["traditional_score_scaled"] = minmax(con["rank_score"])
    con["ml_score_scaled"] = minmax(con["combined_ml_biomarker_score"])
    con["rf_score_scaled"] = minmax(con["random_forest_importance"])
    con["consensus_score"] = con[
        ["traditional_score_scaled", "ml_score_scaled", "rf_score_scaled"]
    ].mean(axis=1)
    con["mean_rank"] = con[["traditional_rank", "ml_rank", "rf_rank"]].mean(axis=1)
    return con.sort_values(["consensus_score", "mean_rank"], ascending=[False, True]).reset_index(drop=True)


def save_pseudobulk(
    mean_expr: np.ndarray,
    meta: pd.DataFrame,
    genes: np.ndarray,
    path: Path,
) -> None:
    expr = pd.DataFrame(mean_expr, columns=genes)
    pd.concat([meta.reset_index(drop=True), expr], axis=1).to_csv(path, index=False)


def plot_bar(df: pd.DataFrame, value_col: str, title: str, path: Path, top_n: int, color: str) -> None:
    top = df.sort_values(value_col, ascending=False).head(top_n).iloc[::-1]
    if top.empty:
        return
    fig, ax = plt.subplots(figsize=(7.5, max(5.0, len(top) * 0.28)), dpi=180)
    ax.barh(top["gene"], top[value_col], color=color)
    ax.set_xlabel(value_col.replace("_", " "))
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_volcano(de: pd.DataFrame, title: str, path: Path) -> None:
    df = de.copy()
    df["neg_log10_fdr"] = -np.log10(np.clip(df["adj_p_value"], 1e-300, 1.0))
    sig_up = (df["adj_p_value"] < 0.05) & (df["diff_mean_log_expr"] > 0) & df["passes_expression_filter"]
    sig_down = (df["adj_p_value"] < 0.05) & (df["diff_mean_log_expr"] < 0) & df["passes_expression_filter"]
    colors = np.where(sig_up, "#b24a3b", np.where(sig_down, "#2f6f8f", "#b8b8b8"))
    fig, ax = plt.subplots(figsize=(7, 5), dpi=180)
    ax.scatter(df["diff_mean_log_expr"], df["neg_log10_fdr"], s=5, c=colors, alpha=0.65, linewidths=0)
    ax.axhline(-np.log10(0.05), color="#555555", linewidth=0.8, linestyle=":")
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_xlabel("Mean log-expression difference (CKD - Ref)")
    ax.set_ylabel("-log10 FDR")
    ax.set_title(title)
    for _, row in pd.concat([df.loc[sig_up].head(8), df.loc[sig_down].head(8)]).iterrows():
        ax.text(row["diff_mean_log_expr"], row["neg_log10_fdr"], row["gene"], fontsize=6)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_consensus_scores(con: pd.DataFrame, cell_class: str, path: Path, top_n: int) -> None:
    top = con.head(top_n).iloc[::-1]
    if top.empty:
        return
    y = np.arange(len(top))
    height = 0.24
    fig, ax = plt.subplots(figsize=(9, max(5.0, len(top) * 0.30)), dpi=180)
    ax.barh(y - height, top["traditional_score_scaled"], height=height, label="Traditional DE", color="#2f6f8f")
    ax.barh(y, top["ml_score_scaled"], height=height, label="RF+KNN ML", color="#b24a3b")
    ax.barh(y + height, top["rf_score_scaled"], height=height, label="Random Forest", color="#6f8f3a")
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"])
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Scaled score within consensus genes")
    ax.set_title(f"{cell_class}: consensus biomarker score comparison")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_consensus_heatmap(con: pd.DataFrame, cell_class: str, path: Path, top_n: int) -> None:
    top = con.head(top_n).set_index("gene")[
        ["traditional_score_scaled", "ml_score_scaled", "rf_score_scaled"]
    ].rename(
        columns={
            "traditional_score_scaled": "Traditional DE",
            "ml_score_scaled": "RF+KNN ML",
            "rf_score_scaled": "Random Forest",
        }
    )
    top = top.iloc[::-1]
    if top.empty:
        return
    fig, ax = plt.subplots(figsize=(6.5, max(5.0, len(top) * 0.30)), dpi=180)
    image = ax.imshow(top.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(np.arange(top.shape[1]))
    ax.set_xticklabels(top.columns, rotation=25, ha="right")
    ax.set_yticks(np.arange(top.shape[0]))
    ax.set_yticklabels(top.index)
    ax.set_title(f"{cell_class}: consensus evidence heatmap")
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Scaled score")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_cv_metrics(cv_df: pd.DataFrame, title: str, path: Path) -> None:
    metric = "balanced_accuracy"
    sub = cv_df[cv_df["metric"] == metric].copy()
    fig, ax = plt.subplots(figsize=(5.5, 4), dpi=180)
    ax.bar(sub["model"], sub["mean"], yerr=sub["std"], color=["#2f6f8f", "#b24a3b"], capsize=4)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Cross-validated balanced accuracy")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_overlap_counts(overlap_counts: pd.DataFrame, cell_class: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4), dpi=180)
    ax.bar(overlap_counts["top_n"].astype(str), overlap_counts["triple_overlap_genes"], color="#2f6f8f")
    ax.set_xlabel("Top N genes from each method")
    ax.set_ylabel("Triple-overlap genes")
    ax.set_title(f"{cell_class}: DE + ML + RF overlap")
    for i, value in enumerate(overlap_counts["triple_overlap_genes"]):
        ax.text(i, value + max(overlap_counts["triple_overlap_genes"].max() * 0.02, 0.3), str(value), ha="center")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def overlap_counts(de_up: pd.DataFrame, ml_up: pd.DataFrame, rf: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for n in [25, 50, 100, 200, 500]:
        overlap = (
            set(de_up.head(min(n, len(de_up)))["gene"])
            & set(ml_up.head(min(n, len(ml_up)))["gene"])
            & set(rf.head(min(n, len(rf)))["gene"])
        )
        rows.append({"top_n": n, "triple_overlap_genes": len(overlap)})
    rows.append(
        {
            "top_n": "all",
            "triple_overlap_genes": len(set(de_up["gene"]) & set(ml_up["gene"]) & set(rf["gene"])),
        }
    )
    return pd.DataFrame(rows)


def plot_summary_counts(summary: pd.DataFrame, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=180)
    ax.bar(summary["cell_class"], summary["n_consensus_genes"], color="#2f6f8f")
    ax.set_ylabel("Consensus gene count")
    ax.set_title("Cell-type consensus biomarker counts")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(figdir / "cell_type_consensus_counts.png")
    plt.close(fig)


def plot_summary_heatmap(all_top: pd.DataFrame, figdir: Path) -> None:
    if all_top.empty:
        return
    table = all_top.pivot_table(
        index="gene",
        columns="cell_class",
        values="consensus_score",
        aggfunc="max",
        fill_value=0,
    )
    table["max_score"] = table.max(axis=1)
    table = table.sort_values("max_score", ascending=False).drop(columns="max_score").head(40)
    fig, ax = plt.subplots(figsize=(7.5, max(6, table.shape[0] * 0.24)), dpi=180)
    image = ax.imshow(table.values, aspect="auto", cmap="magma", vmin=0, vmax=max(table.values.max(), 1e-9))
    ax.set_xticks(np.arange(table.shape[1]))
    ax.set_xticklabels(table.columns, rotation=25, ha="right")
    ax.set_yticks(np.arange(table.shape[0]))
    ax.set_yticklabels(table.index)
    ax.set_title("Top cell-type consensus biomarkers")
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Consensus score")
    fig.tight_layout()
    fig.savefig(figdir / "cell_type_top_consensus_heatmap.png")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    datadir = Path(args.datadir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    datadir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input, backed="r")
    obs = adata.obs.copy()
    for col in [args.condition_col, args.patient_col, args.class_col]:
        obs[col] = obs[col].astype(str)
    matrix, genes = get_expression_matrix_and_genes(adata)
    print(f"Using expression matrix with {len(genes):,} genes")

    summaries: list[dict[str, object]] = []
    top_consensus_rows: list[pd.DataFrame] = []

    for cell_class in DEFAULT_CLASSES:
        slug = slugify(cell_class)
        sums, detected, meta = aggregate_class_patient(
            matrix,
            obs,
            genes,
            cell_class,
            args.patient_col,
            args.condition_col,
            args.class_col,
            args.case,
            args.control,
            args.chunk_size,
        )
        keep_patient = meta["n_cells"].to_numpy() >= args.min_cells_per_patient
        sums = sums[keep_patient]
        detected = detected[keep_patient]
        meta = meta.loc[keep_patient].reset_index(drop=True)
        mean_expr = sums / np.maximum(meta["n_cells"].to_numpy()[:, None], 1)

        save_pseudobulk(mean_expr, meta, genes, datadir / f"{slug}_patient_pseudobulk.csv")
        meta.to_csv(datadir / f"{slug}_patient_metadata.csv", index=False)

        de = compare_de(
            mean_expr,
            detected,
            meta,
            genes,
            args.condition_col,
            args.case,
            args.control,
            args.min_pct,
            args.min_mean,
        )
        de.insert(0, "cell_class", cell_class)
        de.to_csv(outdir / f"{slug}_traditional_de.csv", index=False)
        de_up = de[(de["diff_mean_log_expr"] > 0) & de["passes_expression_filter"]].copy()
        de_sig_up = de_up[de_up["adj_p_value"] < 0.05].copy()
        de_down = de[(de["diff_mean_log_expr"] < 0) & de["passes_expression_filter"]].copy()
        de_up.to_csv(outdir / f"{slug}_traditional_up_biomarkers.csv", index=False)
        de_sig_up.to_csv(outdir / f"{slug}_traditional_significant_up_biomarkers.csv", index=False)
        de_down.to_csv(outdir / f"{slug}_traditional_down_genes.csv", index=False)

        print(f"{cell_class}: running Random Forest and KNN")
        rf, knn_perm, knn_single, cv_metrics, final_fit = fit_ml(
            mean_expr,
            meta,
            genes,
            args.condition_col,
            args.case,
            args.control,
            args.top_n_filter,
            args.rf_trees,
            args.permutation_repeats,
            args.random_state,
        )
        rf.to_csv(outdir / f"{slug}_random_forest_importance.csv", index=False)
        knn_perm.to_csv(outdir / f"{slug}_knn_permutation_importance.csv", index=False)
        knn_single.to_csv(outdir / f"{slug}_knn_single_gene_importance.csv", index=False)
        cv_metrics.insert(0, "cell_class", cell_class)
        cv_metrics.to_csv(outdir / f"{slug}_model_cv_metrics.csv", index=False)
        final_fit.insert(0, "cell_class", cell_class)
        final_fit.to_csv(outdir / f"{slug}_final_fit_diagnostics.csv", index=False)

        ml = combine_ml(de, rf, knn_perm, knn_single)
        ml.insert(0, "cell_class", cell_class)
        ml.to_csv(outdir / f"{slug}_ml_combined_biomarkers.csv", index=False)
        ml_up = ml[ml["diff_mean_log_expr"] > 0].copy()
        ml_up.to_csv(outdir / f"{slug}_ml_upregulated_biomarkers.csv", index=False)

        con = consensus_table(de_sig_up, ml_up, rf)
        con.insert(0, "cell_class", cell_class)
        con.to_csv(outdir / f"{slug}_consensus_traditional_ml_rf.csv", index=False)
        oc = overlap_counts(de_sig_up, ml_up, rf)
        oc.insert(0, "cell_class", cell_class)
        oc.to_csv(outdir / f"{slug}_consensus_overlap_counts.csv", index=False)

        plot_volcano(de, f"{cell_class}: CKD vs Ref traditional DE", figdir / f"{slug}_traditional_volcano.png")
        plot_bar(de_up, "rank_score", f"{cell_class}: top CKD-up traditional biomarkers", figdir / f"{slug}_traditional_top20_up.png", args.plot_top_n, "#2f6f8f")
        plot_bar(rf, "random_forest_importance", f"{cell_class}: Random Forest importance", figdir / f"{slug}_random_forest_top20.png", args.plot_top_n, "#6f8f3a")
        plot_bar(knn_single, "knn_single_gene_importance_score", f"{cell_class}: KNN single-gene importance", figdir / f"{slug}_knn_single_gene_top20.png", args.plot_top_n, "#8060a8")
        plot_bar(ml_up, "combined_ml_biomarker_score", f"{cell_class}: top CKD-up ML biomarkers", figdir / f"{slug}_ml_top20_up.png", args.plot_top_n, "#b24a3b")
        plot_consensus_scores(con, cell_class, figdir / f"{slug}_consensus_score_comparison.png", args.plot_top_n)
        plot_consensus_heatmap(con, cell_class, figdir / f"{slug}_consensus_score_heatmap.png", args.plot_top_n)
        plot_overlap_counts(oc, cell_class, figdir / f"{slug}_consensus_overlap_counts.png")
        plot_cv_metrics(cv_metrics, f"{cell_class}: model CV performance", figdir / f"{slug}_model_cv_balanced_accuracy.png")

        top_part = con.head(10).copy()
        if not top_part.empty:
            top_consensus_rows.append(top_part)

        summaries.append(
            {
                "cell_class": cell_class,
                "n_patients": meta.shape[0],
                "n_case_patients": int((meta[args.condition_col] == args.case).sum()),
                "n_control_patients": int((meta[args.condition_col] == args.control).sum()),
                "n_cells": int(meta["n_cells"].sum()),
                "n_traditional_up": de_up.shape[0],
                "n_traditional_significant_up": de_sig_up.shape[0],
                "n_ml_up": ml_up.shape[0],
                "n_consensus_genes": con.shape[0],
                "top_consensus_genes": ";".join(con["gene"].head(10).tolist()),
            }
        )

        print(f"{cell_class}: top consensus genes")
        print(con[["gene", "consensus_score", "traditional_rank", "ml_rank", "rf_rank"]].head(10).to_string(index=False))

    summary = pd.DataFrame(summaries)
    summary.to_csv(outdir / "cell_type_biomarker_summary.csv", index=False)
    all_top = pd.concat(top_consensus_rows, axis=0, ignore_index=True) if top_consensus_rows else pd.DataFrame()
    if not all_top.empty:
        all_top.to_csv(outdir / "cell_type_top_consensus_genes.csv", index=False)
    plot_summary_counts(summary, figdir)
    plot_summary_heatmap(all_top, figdir)
    adata.file.close()

    print("\nDone. Summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
