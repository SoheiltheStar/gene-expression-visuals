from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from gene_expression_visuals import analysis, data_loader, ml, visualization


st.set_page_config(page_title="Gene Expression Insights", layout="wide")


@st.cache_data(show_spinner=False)
def load_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    differential = data_loader.load_differential_expression()
    expression = data_loader.load_expression_long()
    metadata = data_loader.load_metadata()
    return differential, expression, metadata


de_table, expression_long, metadata = load_tables()


st.title("Airway RNA-seq insight explorer")
st.write(
    "Interact with a compact story derived from the GEO airway dataset (GSE5278)."
)

with st.sidebar:
    st.header("Controls")
    padj_threshold = st.slider(
        "Adjusted p-value cutoff",
        min_value=0.0001,
        max_value=0.25,
        value=0.05,
        step=0.0001,
        format="%.4f",
    )
    log2_fc_threshold = st.slider(
        "Absolute log2 fold-change",
        min_value=0.5,
        max_value=4.0,
        value=1.0,
        step=0.1,
    )
    heatmap_genes = st.slider(
        "Genes in heatmap",
        min_value=10,
        max_value=60,
        value=30,
        step=5,
    )
    classifier_genes = st.slider(
        "Top genes for classifier",
        min_value=8,
        max_value=60,
        value=20,
        step=4,
    )


significant = analysis.filter_differential_expression(
    de_table, padj_threshold, log2_fc_threshold
)
st.metric("Significant genes", f"{len(significant):,}")

volcano = visualization.make_volcano_plot(
    de_table,
    padj_threshold=padj_threshold,
    log2_fc_threshold=log2_fc_threshold,
)
st.subheader("Differential expression volcano plot")
st.plotly_chart(volcano, use_container_width=True)
st.dataframe(significant.head(25))


selected_heatmap_genes = analysis.top_genes_by_significance(de_table, heatmap_genes)
heatmap_matrix = analysis.prepare_heatmap_data(expression_long, selected_heatmap_genes)
st.subheader("Clustered heatmap of top genes")
st.plotly_chart(visualization.make_heatmap(heatmap_matrix), use_container_width=True)


pca_scores = analysis.compute_pca_scores(expression_long, metadata)
st.subheader("Sample separation via PCA")
st.plotly_chart(visualization.make_pca_scatter(pca_scores), use_container_width=True)


classifier_geneset = analysis.top_genes_by_significance(de_table, classifier_genes)
metrics, importances, roc_points = ml.run_logistic_classifier(
    expression_long,
    metadata,
    classifier_geneset,
)
col1, col2 = st.columns([2, 3])
with col1:
    st.subheader("Classifier metrics")
    st.dataframe(metrics.style.format({"value": "{:.3f}"}))
    st.subheader("Top feature importances")
    st.dataframe(importances.head(15).style.format({"importance": "{:.3f}"}))
with col2:
    st.subheader("Classifier diagnostics")
    st.plotly_chart(visualization.make_feature_importance_bar(importances), use_container_width=True)
    st.plotly_chart(visualization.make_roc_curve(roc_points), use_container_width=True)


st.subheader("Peek at expression table")
st.dataframe(expression_long.head(50))

