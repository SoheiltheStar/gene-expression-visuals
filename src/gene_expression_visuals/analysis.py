"""Analysis helpers for interactive visualisations."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def filter_differential_expression(
    de_table: pd.DataFrame,
    padj_max: float,
    log2_fc_min: float,
) -> pd.DataFrame:
    mask = (de_table["padj"] <= padj_max) & (
        de_table["log2_fold_change"].abs() >= log2_fc_min
    )
    return de_table.loc[mask].copy()


def top_genes_by_significance(de_table: pd.DataFrame, limit: int) -> pd.Series:
    sorted_df = de_table.sort_values(["padj", "log2_fold_change"], ascending=[True, False])
    return sorted_df.head(limit)["gene"]


def pivot_expression(expression_long: pd.DataFrame) -> pd.DataFrame:
    matrix = expression_long.pivot_table(
        index="sample_id",
        columns="gene",
        values="log2_expression",
    )
    return matrix.sort_index(axis=0).sort_index(axis=1)


def prepare_heatmap_data(
    expression_long: pd.DataFrame,
    genes: Iterable[str],
) -> pd.DataFrame:
    subset = expression_long[expression_long["gene"].isin(list(genes))]
    matrix = subset.pivot_table(
        index="gene",
        columns="sample_id",
        values="log2_expression",
    )
    return matrix


def compute_pca_scores(
    expression_long: pd.DataFrame,
    metadata: pd.DataFrame,
    n_components: int = 2,
) -> pd.DataFrame:
    matrix = pivot_expression(expression_long)
    matched_metadata = metadata.set_index("sample_id").loc[matrix.index]
    pca = PCA(n_components=n_components, random_state=42)
    scores = pca.fit_transform(matrix.values)
    columns = [f"PC{i+1}" for i in range(n_components)]
    df = pd.DataFrame(scores, columns=columns, index=matrix.index).reset_index()
    df = df.merge(matched_metadata.reset_index(), on="sample_id", how="left")
    variance = np.round(pca.explained_variance_ratio_, 4)
    for idx, value in enumerate(variance, start=1):
        df[f"explained_variance_pc{idx}"] = value
    return df

