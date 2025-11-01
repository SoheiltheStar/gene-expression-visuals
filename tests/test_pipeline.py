from pathlib import Path

from gene_expression_visuals import analysis, data_loader, ml


DATA_DIR = Path("data/processed")


def test_processed_tables_exist() -> None:
    de = data_loader.load_differential_expression(DATA_DIR)
    expression = data_loader.load_expression_long(DATA_DIR)
    metadata = data_loader.load_metadata(DATA_DIR)
    assert not de.empty
    assert not expression.empty
    assert not metadata.empty


def test_analysis_helpers() -> None:
    de = data_loader.load_differential_expression(DATA_DIR)
    expression = data_loader.load_expression_long(DATA_DIR)
    metadata = data_loader.load_metadata(DATA_DIR)
    filtered = analysis.filter_differential_expression(de, padj_max=0.1, log2_fc_min=0.5)
    assert filtered["padj"].max() <= 0.1
    top_genes = analysis.top_genes_by_significance(de, limit=20)
    heatmap_matrix = analysis.prepare_heatmap_data(expression, top_genes)
    pca_scores = analysis.compute_pca_scores(expression, metadata)
    assert len(top_genes) == 20
    assert heatmap_matrix.shape[0] <= 20
    assert {"PC1", "PC2"}.issubset(pca_scores.columns)


def test_classifier_pipeline() -> None:
    de = data_loader.load_differential_expression(DATA_DIR)
    expression = data_loader.load_expression_long(DATA_DIR)
    metadata = data_loader.load_metadata(DATA_DIR)
    genes = analysis.top_genes_by_significance(de, limit=25)
    metrics, importances, roc_points = ml.run_logistic_classifier(expression, metadata, genes)
    assert metrics.loc[metrics["metric"] == "accuracy", "value"].iloc[0] > 0.5
    assert not importances.empty
    assert not roc_points.empty

