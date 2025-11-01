"""Lightweight machine-learning helpers for classifier experiments."""

from __future__ import annotations

from typing import Iterable, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .analysis import pivot_expression


def encode_condition(series: pd.Series) -> pd.Series:
    values = series.str.lower().map({"treated": 1, "control": 0, "untreated": 0})
    if values.isna().any():
        raise ValueError("Condition column must contain only treated/control labels")
    return values


def prepare_feature_matrix(
    expression_long: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: Iterable[str],
) -> Tuple[pd.DataFrame, pd.Series]:
    matrix = pivot_expression(expression_long)
    subset = matrix[list(genes)]
    aligned_metadata = metadata.set_index("sample_id").loc[subset.index]
    labels = encode_condition(aligned_metadata["condition"])
    return subset, labels


def run_logistic_classifier(
    expression_long: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: Iterable[str],
    test_size: float = 0.25,
    random_state: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    features, labels = prepare_feature_matrix(expression_long, metadata, genes)
    train_x, test_x, train_y, test_y = train_test_split(
        features,
        labels,
        test_size=test_size,
        random_state=random_state,
        stratify=labels,
    )
    pipeline = Pipeline(
        steps=[
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    penalty="l2",
                    solver="liblinear",
                    random_state=random_state,
                ),
            ),
        ]
    )
    pipeline.fit(train_x, train_y)
    probs = pipeline.predict_proba(test_x)[:, 1]
    preds = pipeline.predict(test_x)
    metrics = pd.DataFrame(
        {
            "metric": ["accuracy", "roc_auc"],
            "value": [accuracy_score(test_y, preds), roc_auc_score(test_y, probs)],
        }
    )
    coef = pipeline.named_steps["clf"].coef_[0]
    importances = (
        pd.DataFrame({"gene": list(genes), "importance": np.abs(coef)})
        .sort_values("importance", ascending=False)
        .reset_index(drop=True)
    )
    fpr, tpr, thresholds = roc_curve(test_y, probs)
    roc_points = pd.DataFrame({"fpr": fpr, "tpr": tpr, "threshold": thresholds})
    return metrics, importances, roc_points

