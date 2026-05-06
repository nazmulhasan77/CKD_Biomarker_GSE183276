"""
Build patient-level pseudobulk features for ML biomarker discovery.

Input:
  ../GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad

Output:
  Using ML/data/patient_pseudobulk_logexpr.csv
  Using ML/data/patient_metadata.csv

Why patient-level?
  The patient is the biological replicate. Training ML on individual cells would
  leak patient identity and overstate model performance.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


DEFAULT_INPUT = "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5ad"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create patient pseudobulk matrix.")
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Input .h5ad file")
    parser.add_argument("--outdir", default="Using ML/data", help="Output data folder")
    parser.add_argument("--condition-col", default="condition.l1")
    parser.add_argument("--patient-col", default="patient")
    parser.add_argument("--case", default="CKD")
    parser.add_argument("--control", default="Ref")
    parser.add_argument("--chunk-size", type=int, default=2000)
    return parser.parse_args()


def make_unique_gene_names(genes: np.ndarray) -> np.ndarray:
    seen: dict[str, int] = {}
    out: list[str] = []
    for gene in genes.astype(str):
        count = seen.get(gene, 0)
        out.append(gene if count == 0 else f"{gene}__dup{count}")
        seen[gene] = count + 1
    return np.array(out, dtype=object)


def get_raw_matrix_and_genes(adata: ad.AnnData) -> tuple[object, np.ndarray]:
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


def aggregate_patient_means(
    matrix: object,
    obs: pd.DataFrame,
    genes: np.ndarray,
    patient_col: str,
    condition_col: str,
    case: str,
    control: str,
    chunk_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    keep = obs[condition_col].isin([case, control]).to_numpy()
    selected = np.flatnonzero(keep)
    if selected.size == 0:
        raise ValueError("No CKD/Ref cells found.")

    group_obs = obs.iloc[selected][[patient_col, condition_col]].copy()
    group_obs[patient_col] = group_obs[patient_col].astype(str)
    group_obs[condition_col] = group_obs[condition_col].astype(str)
    keys = (group_obs[patient_col] + "||" + group_obs[condition_col]).to_numpy()

    unique_keys, inverse = np.unique(keys, return_inverse=True)
    n_groups = unique_keys.size
    n_genes = len(genes)

    sums = np.zeros((n_groups, n_genes), dtype=np.float64)
    n_cells = np.zeros(n_groups, dtype=np.int64)

    print(f"Aggregating {selected.size:,} cells into {n_groups} patient profiles.")
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
            n_cells[code] += local.size

        if start == 0 or stop == selected.size or (stop // chunk_size) % 20 == 0:
            print(f"  processed {stop:,}/{selected.size:,} cells")

    means = sums / np.maximum(n_cells[:, None], 1)
    metadata = pd.DataFrame(
        [key.split("||") for key in unique_keys], columns=[patient_col, condition_col]
    )
    metadata["n_cells"] = n_cells

    expr = pd.DataFrame(means, columns=genes)
    pseudobulk = pd.concat([metadata, expr], axis=1)
    pseudobulk = pseudobulk.sort_values([condition_col, patient_col]).reset_index(drop=True)
    metadata = pseudobulk[[patient_col, condition_col, "n_cells"]].copy()
    return pseudobulk, metadata


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.input}")
    adata = ad.read_h5ad(args.input, backed="r")
    obs = adata.obs.copy()
    matrix, genes = get_raw_matrix_and_genes(adata)
    print(f"Using expression matrix with {len(genes):,} genes.")

    pseudobulk, metadata = aggregate_patient_means(
        matrix=matrix,
        obs=obs,
        genes=genes,
        patient_col=args.patient_col,
        condition_col=args.condition_col,
        case=args.case,
        control=args.control,
        chunk_size=args.chunk_size,
    )

    pseudobulk_path = outdir / "patient_pseudobulk_logexpr.csv"
    metadata_path = outdir / "patient_metadata.csv"
    pseudobulk.to_csv(pseudobulk_path, index=False)
    metadata.to_csv(metadata_path, index=False)
    adata.file.close()

    print(f"Saved: {pseudobulk_path}")
    print(f"Saved: {metadata_path}")
    print(metadata[args.condition_col].value_counts().to_string())


if __name__ == "__main__":
    main()
