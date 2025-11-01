"""Plotly figure factories for the gene expression dashboard."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def make_volcano_plot(
    de_table: pd.DataFrame,
    padj_threshold: float,
    log2_fc_threshold: float,
) -> go.Figure:
    df = de_table.copy()
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["significant"] = (
        (df["padj"] <= padj_threshold)
        & (df["log2_fold_change"].abs() >= log2_fc_threshold)
    )
    fig = px.scatter(
        df,
        x="log2_fold_change",
        y="neg_log10_padj",
        color="significant",
        hover_data={"gene": True, "padj": ":.2e", "log2_fold_change": ":.2f"},
        labels={
            "log2_fold_change": "log2 fold change",
            "neg_log10_padj": "-log10 adjusted p-value",
        },
        color_discrete_map={True: "crimson", False: "steelblue"},
    )
    fig.update_layout(title="Differential expression volcano plot", legend_title="Significant")
    fig.add_vline(x=log2_fc_threshold, line_dash="dash", line_color="gray")
    fig.add_vline(x=-log2_fc_threshold, line_dash="dash", line_color="gray")
    fig.add_hline(y=-np.log10(padj_threshold), line_dash="dash", line_color="gray")
    return fig


def make_heatmap(matrix: pd.DataFrame) -> go.Figure:
    fig = px.imshow(
        matrix,
        color_continuous_scale="Viridis",
        aspect="auto",
        labels={"x": "Sample", "y": "Gene", "color": "log2 expression"},
    )
    fig.update_layout(title="Top differentially expressed genes heatmap")
    return fig


def make_pca_scatter(pca_scores: pd.DataFrame) -> go.Figure:
    explained = (
        pca_scores[["explained_variance_pc1", "explained_variance_pc2"]]
        .iloc[0]
        .to_dict()
    )
    axis_labels = {
        "PC1": f"PC1 ({explained['explained_variance_pc1']*100:.1f}% var)",
        "PC2": f"PC2 ({explained['explained_variance_pc2']*100:.1f}% var)",
    }
    fig = px.scatter(
        pca_scores,
        x="PC1",
        y="PC2",
        color="condition",
        hover_data={"sample_id": True, "title": True},
        labels=axis_labels,
    )
    fig.update_layout(title="Principal component view of samples")
    return fig


def make_feature_importance_bar(importances: pd.DataFrame) -> go.Figure:
    fig = px.bar(
        importances.head(15),
        x="importance",
        y="gene",
        orientation="h",
        labels={"importance": "Coefficient magnitude", "gene": ""},
    )
    fig.update_layout(title="Classifier feature importances", yaxis=dict(autorange="reversed"))
    return fig


def make_roc_curve(roc_points: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=roc_points["fpr"],
            y=roc_points["tpr"],
            mode="lines",
            name="Model",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode="lines",
            name="Chance",
            line=dict(dash="dash", color="gray"),
        )
    )
    fig.update_layout(
        title="ROC curve",
        xaxis_title="False positive rate",
        yaxis_title="True positive rate",
        yaxis=dict(scaleanchor="x", scaleratio=1),
    )
    return fig

