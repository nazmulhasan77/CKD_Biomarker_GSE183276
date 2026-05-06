"""
Random Forest and KNN biomarker discovery from patient-level pseudobulk data.

Outputs:
  Using ML/results/random_forest_gene_importance.csv
  Using ML/results/knn_permutation_importance.csv
  Using ML/results/knn_single_gene_cv_importance.csv
  Using ML/results/ml_combined_candidate_biomarkers.csv
  Using ML/results/ml_upregulated_candidate_biomarkers.csv
  Using ML/results/ml_downregulated_reference_loss_genes.csv
  Using ML/results/model_cv_metrics.csv
  Using ML/figures/*.png

KNN has no native feature-importance value, so this script uses permutation
importance after fitting the best KNN model.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_validate
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

warnings.filterwarnings("ignore", message="Features .* are constant", category=UserWarning)
warnings.filterwarnings("ignore", message="invalid value encountered in divide", category=RuntimeWarning)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find CKD candidate biomarkers with Random Forest and KNN."
    )
    parser.add_argument(
        "--input",
        default="Using ML/data/patient_pseudobulk_logexpr.csv",
        help="Patient pseudobulk CSV from 01_build_patient_pseudobulk.py",
    )
    parser.add_argument("--outdir", default="Using ML/results")
    parser.add_argument("--figdir", default="Using ML/figures")
    parser.add_argument("--condition-col", default="condition.l1")
    parser.add_argument("--patient-col", default="patient")
    parser.add_argument("--case", default="CKD")
    parser.add_argument("--control", default="Ref")
    parser.add_argument(
        "--top-n-filter",
        type=int,
        default=500,
        help="Top univariate genes retained before fitting final ML models.",
    )
    parser.add_argument("--rf-trees", type=int, default=1000)
    parser.add_argument("--permutation-repeats", type=int, default=50)
    parser.add_argument("--random-state", type=int, default=7)
    return parser.parse_args()


def safe_f_classif(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Warning-free ANOVA F score for SelectKBest."""
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
        if Xc.shape[0] == 0:
            continue
        class_mean = Xc.mean(axis=0)
        ss_between += Xc.shape[0] * (class_mean - overall_mean) ** 2
        ss_within += ((Xc - class_mean) ** 2).sum(axis=0)

    df_between = max(n_classes - 1, 1)
    df_within = max(n_samples - n_classes, 1)
    ms_between = ss_between / df_between
    ms_within = ss_within / df_within
    scores = np.divide(
        ms_between,
        ms_within,
        out=np.zeros_like(ms_between),
        where=ms_within > 0,
    )
    p_values = stats.f.sf(scores, df_between, df_within)
    scores = np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)
    p_values = np.nan_to_num(p_values, nan=1.0, posinf=1.0, neginf=1.0)
    return scores, p_values


def load_data(
    path: Path, condition_col: str, patient_col: str, case: str, control: str
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(path)
    required = {patient_col, condition_col, "n_cells"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    df = df[df[condition_col].isin([case, control])].copy()
    df = df.sort_values([condition_col, patient_col]).reset_index(drop=True)
    gene_cols = [
        col for col in df.columns if col not in {patient_col, condition_col, "n_cells"}
    ]
    X = df[gene_cols].to_numpy(dtype=np.float64)
    y = (df[condition_col].to_numpy() == case).astype(int)
    genes = np.array(gene_cols, dtype=object)

    variance = X.var(axis=0)
    keep = variance > 0
    X = X[:, keep]
    genes = genes[keep]
    return df[[patient_col, condition_col, "n_cells"]], X, y, genes


def patient_level_stats(
    X: np.ndarray, y: np.ndarray, genes: np.ndarray, case: str, control: str
) -> pd.DataFrame:
    case_X = X[y == 1]
    control_X = X[y == 0]
    case_mean = case_X.mean(axis=0)
    control_mean = control_X.mean(axis=0)
    diff = case_mean - control_mean

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        t_stat, p_value = stats.ttest_ind(
            case_X, control_X, axis=0, equal_var=False, nan_policy="omit"
        )
    p_value = np.nan_to_num(p_value, nan=1.0, posinf=1.0, neginf=1.0)
    adj_p = multipletests(p_value, method="fdr_bh")[1]

    n_case = case_X.shape[0]
    n_control = control_X.shape[0]
    case_sd = case_X.std(axis=0, ddof=1)
    control_sd = control_X.std(axis=0, ddof=1)
    pooled_sd = np.sqrt(
        ((n_case - 1) * case_sd**2 + (n_control - 1) * control_sd**2)
        / max(n_case + n_control - 2, 1)
    )
    cohen_d = np.divide(
        diff,
        pooled_sd,
        out=np.zeros_like(diff),
        where=pooled_sd > 0,
    )

    f_score, f_p = safe_f_classif(X, y)
    stat_score = -np.log10(np.clip(adj_p, 1e-300, 1.0)) * np.abs(cohen_d) * np.abs(diff)

    return pd.DataFrame(
        {
            "gene": genes,
            f"mean_expr_{case}": case_mean,
            f"mean_expr_{control}": control_mean,
            "diff_mean_log_expr": diff,
            "direction": np.where(diff >= 0, f"up_in_{case}", f"down_in_{case}"),
            "cohen_d": cohen_d,
            "t_stat": t_stat,
            "p_value": p_value,
            "adj_p_value": adj_p,
            "f_classif_score": f_score,
            "f_classif_p_value": f_p,
            "statistical_rank_score": stat_score,
        }
    )


def make_cv(y: np.ndarray, random_state: int) -> StratifiedKFold:
    class_counts = np.bincount(y)
    n_splits = int(min(5, class_counts.min()))
    if n_splits < 3:
        raise ValueError("Need at least 3 patients per class for CV.")
    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)


def summarize_cross_validation(
    X: np.ndarray,
    y: np.ndarray,
    top_n_filter: int,
    rf_trees: int,
    random_state: int,
) -> pd.DataFrame:
    k = min(top_n_filter, X.shape[1])
    outer_cv = make_cv(y, random_state)
    inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=random_state)

    rf_pipe = Pipeline(
        steps=[
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
    )

    knn_pipe = Pipeline(
        steps=[
            ("select", SelectKBest(score_func=safe_f_classif, k=k)),
            ("scale", StandardScaler()),
            ("model", KNeighborsClassifier()),
        ]
    )
    knn_grid = GridSearchCV(
        estimator=knn_pipe,
        param_grid={
            "model__n_neighbors": [3, 5, 7, 9, 11],
            "model__weights": ["uniform", "distance"],
        },
        scoring="balanced_accuracy",
        cv=inner_cv,
        n_jobs=-1,
        refit=True,
    )

    scoring = {
        "accuracy": "accuracy",
        "balanced_accuracy": "balanced_accuracy",
        "f1": "f1",
        "roc_auc": "roc_auc",
    }

    rows: list[dict[str, float | str]] = []
    for name, estimator in [("RandomForest", rf_pipe), ("KNN", knn_grid)]:
        scores = cross_validate(
            estimator,
            X,
            y,
            cv=outer_cv,
            scoring=scoring,
            n_jobs=-1,
            error_score=np.nan,
        )
        for metric in scoring:
            values = scores[f"test_{metric}"]
            rows.append(
                {
                    "model": name,
                    "metric": metric,
                    "mean": float(np.nanmean(values)),
                    "std": float(np.nanstd(values, ddof=1)),
                    "fold_values": ";".join(f"{v:.4f}" for v in values),
                }
            )
    return pd.DataFrame(rows)


def final_feature_selection(
    X: np.ndarray, y: np.ndarray, genes: np.ndarray, top_n_filter: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, SelectKBest]:
    k = min(top_n_filter, X.shape[1])
    selector = SelectKBest(score_func=safe_f_classif, k=k)
    selector.fit(X, y)
    idx = selector.get_support(indices=True)
    return X[:, idx], genes[idx], idx, selector


def fit_random_forest(
    X_selected: np.ndarray,
    y: np.ndarray,
    selected_genes: np.ndarray,
    rf_trees: int,
    random_state: int,
) -> tuple[RandomForestClassifier, pd.DataFrame]:
    rf = RandomForestClassifier(
        n_estimators=rf_trees,
        random_state=random_state,
        class_weight="balanced",
        max_features="sqrt",
        n_jobs=-1,
    )
    rf.fit(X_selected, y)
    importance = pd.DataFrame(
        {
            "gene": selected_genes,
            "random_forest_importance": rf.feature_importances_,
        }
    ).sort_values("random_forest_importance", ascending=False)
    return rf, importance.reset_index(drop=True)


def fit_knn_and_permutation_importance(
    X_selected: np.ndarray,
    y: np.ndarray,
    selected_genes: np.ndarray,
    random_state: int,
    repeats: int,
) -> tuple[Pipeline, pd.DataFrame, pd.DataFrame]:
    cv = make_cv(y, random_state)
    pipe = Pipeline(
        steps=[
            ("scale", StandardScaler()),
            ("model", KNeighborsClassifier()),
        ]
    )
    grid = GridSearchCV(
        pipe,
        param_grid={
            "model__n_neighbors": [3, 5, 7, 9, 11],
            "model__weights": ["uniform", "distance"],
        },
        scoring="balanced_accuracy",
        cv=cv,
        n_jobs=-1,
        refit=True,
    )
    grid.fit(X_selected, y)
    best_model = grid.best_estimator_

    perm = permutation_importance(
        best_model,
        X_selected,
        y,
        scoring="balanced_accuracy",
        n_repeats=repeats,
        random_state=random_state,
        n_jobs=-1,
    )
    importance = pd.DataFrame(
        {
            "gene": selected_genes,
            "knn_permutation_importance_mean": perm.importances_mean,
            "knn_permutation_importance_std": perm.importances_std,
        }
    ).sort_values("knn_permutation_importance_mean", ascending=False)

    grid_results = pd.DataFrame(grid.cv_results_)
    grid_results = grid_results[
        [
            "param_model__n_neighbors",
            "param_model__weights",
            "mean_test_score",
            "std_test_score",
            "rank_test_score",
        ]
    ].sort_values("rank_test_score")
    return best_model, importance.reset_index(drop=True), grid_results


def knn_single_gene_cv_importance(
    X_selected: np.ndarray,
    y: np.ndarray,
    selected_genes: np.ndarray,
    random_state: int,
) -> pd.DataFrame:
    """
    Rank genes by single-gene KNN cross-validated performance.

    This is useful because standard KNN has no coefficient/importance attribute,
    and one-at-a-time permutation can be zero when many genes are redundant.
    """
    cv = make_cv(y, random_state)
    pipe = Pipeline(
        steps=[
            ("scale", StandardScaler()),
            ("model", KNeighborsClassifier(n_neighbors=5, weights="distance")),
        ]
    )
    scoring = {
        "accuracy": "accuracy",
        "balanced_accuracy": "balanced_accuracy",
        "roc_auc": "roc_auc",
    }
    rows: list[dict[str, float | str]] = []
    for idx, gene in enumerate(selected_genes):
        X_gene = X_selected[:, [idx]]
        scores = cross_validate(
            pipe,
            X_gene,
            y,
            cv=cv,
            scoring=scoring,
            n_jobs=1,
            error_score=np.nan,
        )
        mean_balanced = float(np.nanmean(scores["test_balanced_accuracy"]))
        rows.append(
            {
                "gene": gene,
                "knn_single_gene_accuracy_mean": float(
                    np.nanmean(scores["test_accuracy"])
                ),
                "knn_single_gene_balanced_accuracy_mean": mean_balanced,
                "knn_single_gene_roc_auc_mean": float(
                    np.nanmean(scores["test_roc_auc"])
                ),
                "knn_single_gene_importance_score": max(mean_balanced - 0.5, 0.0),
            }
        )
    return pd.DataFrame(rows).sort_values(
        "knn_single_gene_importance_score", ascending=False
    ).reset_index(drop=True)


def scaled_positive(series: pd.Series) -> pd.Series:
    values = series.clip(lower=0).astype(float)
    max_value = values.max()
    if max_value <= 0:
        return values
    return values / max_value


def combine_rankings(
    stats_df: pd.DataFrame,
    rf_importance: pd.DataFrame,
    knn_importance: pd.DataFrame,
    knn_single_gene: pd.DataFrame,
) -> pd.DataFrame:
    out = stats_df.merge(rf_importance, on="gene", how="inner")
    out = out.merge(knn_importance, on="gene", how="inner")
    out = out.merge(knn_single_gene, on="gene", how="inner")
    out["rf_importance_scaled"] = scaled_positive(out["random_forest_importance"])
    out["knn_permutation_importance_scaled"] = scaled_positive(
        out["knn_permutation_importance_mean"]
    )
    out["knn_single_gene_importance_scaled"] = scaled_positive(
        out["knn_single_gene_importance_score"]
    )
    out["statistical_score_scaled"] = scaled_positive(out["statistical_rank_score"])
    out["combined_ml_biomarker_score"] = (
        0.40 * out["rf_importance_scaled"]
        + 0.40 * out["knn_single_gene_importance_scaled"]
        + 0.20 * out["statistical_score_scaled"]
    )
    return out.sort_values("combined_ml_biomarker_score", ascending=False).reset_index(
        drop=True
    )


def save_model_diagnostics(
    rf: RandomForestClassifier,
    knn: Pipeline,
    X_selected: np.ndarray,
    y: np.ndarray,
    outdir: Path,
) -> None:
    rows = []
    for name, model in [("RandomForest_final_fit", rf), ("KNN_final_fit", knn)]:
        pred = model.predict(X_selected)
        if hasattr(model, "predict_proba"):
            prob = model.predict_proba(X_selected)[:, 1]
        else:
            prob = np.full_like(y, fill_value=np.nan, dtype=float)
        tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
        rows.append(
            {
                "model": name,
                "accuracy": accuracy_score(y, pred),
                "balanced_accuracy": balanced_accuracy_score(y, pred),
                "f1": f1_score(y, pred),
                "roc_auc": roc_auc_score(y, prob) if np.isfinite(prob).all() else np.nan,
                "tn": tn,
                "fp": fp,
                "fn": fn,
                "tp": tp,
            }
        )
    pd.DataFrame(rows).to_csv(outdir / "final_fit_diagnostics.csv", index=False)


def plot_top_bar(
    df: pd.DataFrame,
    value_col: str,
    title: str,
    path: Path,
    top_n: int = 20,
    color: str = "#2f6f8f",
) -> None:
    top = df.sort_values(value_col, ascending=False).head(top_n).iloc[::-1]
    if top.empty:
        return
    fig, ax = plt.subplots(figsize=(7.5, 6), dpi=180)
    ax.barh(top["gene"], top[value_col], color=color)
    ax.set_xlabel(value_col.replace("_", " "))
    ax.set_ylabel("")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    metadata, X, y, genes = load_data(
        Path(args.input),
        args.condition_col,
        args.patient_col,
        args.case,
        args.control,
    )
    print(
        f"Loaded {X.shape[0]} patient profiles and {X.shape[1]:,} non-constant genes."
    )
    print(metadata[args.condition_col].value_counts().to_string())

    stats_df = patient_level_stats(X, y, genes, args.case, args.control)
    stats_df.to_csv(outdir / "patient_level_gene_statistics.csv", index=False)

    print("Running leakage-aware cross-validation for RF and KNN.")
    cv_metrics = summarize_cross_validation(
        X=X,
        y=y,
        top_n_filter=args.top_n_filter,
        rf_trees=args.rf_trees,
        random_state=args.random_state,
    )
    cv_metrics.to_csv(outdir / "model_cv_metrics.csv", index=False)

    print(f"Selecting top {args.top_n_filter} genes for final ML importance.")
    X_selected, selected_genes, selected_idx, _ = final_feature_selection(
        X, y, genes, args.top_n_filter
    )
    pd.DataFrame({"gene": selected_genes, "original_gene_index": selected_idx}).to_csv(
        outdir / "selected_genes_for_ml.csv", index=False
    )

    print("Fitting final Random Forest.")
    rf, rf_importance = fit_random_forest(
        X_selected,
        y,
        selected_genes,
        args.rf_trees,
        args.random_state,
    )
    rf_importance.to_csv(outdir / "random_forest_gene_importance.csv", index=False)

    print("Fitting final KNN and computing permutation importance.")
    knn, knn_importance, knn_grid_results = fit_knn_and_permutation_importance(
        X_selected,
        y,
        selected_genes,
        args.random_state,
        args.permutation_repeats,
    )
    knn_importance.to_csv(outdir / "knn_permutation_importance.csv", index=False)
    knn_grid_results.to_csv(outdir / "knn_grid_search_results.csv", index=False)

    print("Running single-gene KNN cross-validation importance.")
    knn_single_gene = knn_single_gene_cv_importance(
        X_selected,
        y,
        selected_genes,
        args.random_state,
    )
    knn_single_gene.to_csv(outdir / "knn_single_gene_cv_importance.csv", index=False)

    combined = combine_rankings(
        stats_df, rf_importance, knn_importance, knn_single_gene
    )
    combined.to_csv(outdir / "ml_combined_candidate_biomarkers.csv", index=False)

    up = combined[combined["diff_mean_log_expr"] > 0].copy()
    down = combined[combined["diff_mean_log_expr"] < 0].copy()
    up.to_csv(outdir / "ml_upregulated_candidate_biomarkers.csv", index=False)
    down.to_csv(outdir / "ml_downregulated_reference_loss_genes.csv", index=False)

    save_model_diagnostics(rf, knn, X_selected, y, outdir)

    plot_top_bar(
        rf_importance,
        "random_forest_importance",
        "Random Forest Gene Importance",
        figdir / "top20_random_forest_importance.png",
        color="#2f6f8f",
    )
    plot_top_bar(
        knn_importance,
        "knn_permutation_importance_mean",
        "KNN Permutation Gene Importance",
        figdir / "top20_knn_permutation_importance.png",
        color="#8060a8",
    )
    plot_top_bar(
        knn_single_gene,
        "knn_single_gene_importance_score",
        "KNN Single-Gene CV Importance",
        figdir / "top20_knn_single_gene_cv_importance.png",
        color="#8060a8",
    )
    plot_top_bar(
        combined,
        "combined_ml_biomarker_score",
        "Combined ML Candidate Biomarkers",
        figdir / "top20_combined_ml_biomarkers.png",
        color="#b24a3b",
    )
    plot_top_bar(
        up,
        "combined_ml_biomarker_score",
        "Combined ML CKD-Up Candidate Biomarkers",
        figdir / "top20_ml_upregulated_biomarkers.png",
        color="#b24a3b",
    )

    print("Top combined ML biomarkers:")
    cols = [
        "gene",
        "combined_ml_biomarker_score",
        "random_forest_importance",
        "knn_single_gene_balanced_accuracy_mean",
        "knn_permutation_importance_mean",
        "diff_mean_log_expr",
        "direction",
        "adj_p_value",
        "cohen_d",
    ]
    print(combined[cols].head(25).to_string(index=False))
    print("\nTop ML CKD-upregulated candidate biomarkers:")
    print(up[cols].head(25).to_string(index=False))


if __name__ == "__main__":
    main()
