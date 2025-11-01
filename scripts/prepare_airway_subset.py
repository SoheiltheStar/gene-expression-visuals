from __future__ import annotations

# pyright: reportMissingImports=false

import gzip
from pathlib import Path
from typing import TYPE_CHECKING, List

import numpy as np
import pandas as pd
import requests

if TYPE_CHECKING:
    from scipy import stats as _stats
    from statsmodels.stats.multitest import multipletests as _multipletests


FPKM_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/suppl/GSE52778_All_Sample_FPKM_Matrix.txt.gz"
RAW_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")


def ensure_directories() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)


def download_file(url: str, destination: Path) -> Path:
    if destination.exists():
        return destination
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    destination.write_bytes(response.content)
    return destination


def load_fpkm_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        df = pd.read_csv(handle, sep=r"\s+", engine="python")
    replicate_columns = [
        col
        for col in df.columns
        if col.startswith("Dex_LL") or col.startswith("Untreated_LL")
    ]
    df["gene"] = df["gene_short_name"].where(df["gene_short_name"] != "-", df["tracking_id"])
    subset = df[["gene"] + replicate_columns]
    subset = subset[subset["gene"].notna()]
    subset = subset.groupby("gene", as_index=True).mean()
    return subset


def build_metadata_from_replicates(replicate_columns: List[str]) -> pd.DataFrame:
    records = []
    for column in replicate_columns:
        parts = column.split("_")
        treatment = parts[0]
        replicate = parts[1] if len(parts) > 1 else "NA"
        condition = "treated" if treatment.lower() == "dex" else "control"
        records.append(
            {
                "sample_id": column,
                "treatment": treatment,
                "replicate": replicate,
                "condition": condition,
            }
        )
    return pd.DataFrame.from_records(records)


def compute_differential_expression(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    control_samples = metadata.loc[metadata["condition"] == "control", "sample_id"]
    treated_samples = metadata.loc[metadata["condition"] == "treated", "sample_id"]
    control = matrix[control_samples]
    treated = matrix[treated_samples]
    log_control = np.log2(control + 1)
    log_treated = np.log2(treated + 1)
    fold_change = log_treated.mean(axis=1) - log_control.mean(axis=1)
    t_stat, p_values = stats.ttest_ind(
        log_treated.values,
        log_control.values,
        axis=1,
        equal_var=False,
    )
    p_values = np.nan_to_num(p_values, nan=1.0)
    padj = multipletests(p_values, method="fdr_bh")[1]
    return (
        pd.DataFrame(
            {
                "gene": matrix.index,
                "log2_fold_change": fold_change.values,
                "p_value": p_values,
                "padj": padj,
            }
        )
        .sort_values("padj")
        .reset_index(drop=True)
    )


def build_expression_long(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    log_matrix = np.log2(matrix + 1)
    tidy = (
        log_matrix.transpose()
        .reset_index()
        .melt(id_vars="index", var_name="gene", value_name="log2_expression")
        .rename(columns={"index": "sample_id"})
    )
    tidy = tidy.merge(metadata, on="sample_id", how="left")
    return tidy


def write_outputs(
    metadata_replicates: pd.DataFrame,
    expression_long: pd.DataFrame,
    de_table: pd.DataFrame,
) -> None:
    metadata_replicates.to_csv(PROCESSED_DIR / "airway_metadata.csv", index=False)
    expression_long.to_csv(PROCESSED_DIR / "airway_expression_long.csv", index=False)
    de_table.to_csv(PROCESSED_DIR / "airway_differential_expression.csv", index=False)


def main() -> None:
    ensure_directories()
    fpkm_path = download_file(FPKM_URL, RAW_DIR / "GSE52778_All_Sample_FPKM_Matrix.txt.gz")
    matrix = load_fpkm_matrix(fpkm_path)
    replicate_columns = list(matrix.columns)
    metadata_replicates = build_metadata_from_replicates(replicate_columns)
    metadata_replicates = metadata_replicates.sort_values("sample_id").reset_index(drop=True)
    expression_long = build_expression_long(matrix, metadata_replicates)
    de_table = compute_differential_expression(matrix, metadata_replicates)
    write_outputs(metadata_replicates, expression_long, de_table)


if __name__ == "__main__":
    main()

