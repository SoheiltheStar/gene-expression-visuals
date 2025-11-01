"""Data loading utilities for processed airway RNA-seq artifacts."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import pandas as pd


DEFAULT_PROCESSED_DIR = Path(os.getenv("DATA_PROCESSED_DIR", "data/processed"))


def resolve_data_dir(data_dir: Optional[str | Path] = None) -> Path:
    path = Path(data_dir) if data_dir else DEFAULT_PROCESSED_DIR
    if not path.exists():
        raise FileNotFoundError(f"Processed data directory not found: {path}")
    return path


def load_differential_expression(data_dir: Optional[str | Path] = None) -> pd.DataFrame:
    base = resolve_data_dir(data_dir)
    table = base / "airway_differential_expression.csv"
    if not table.exists():
        raise FileNotFoundError(f"Missing differential expression table at {table}")
    df = pd.read_csv(table)
    expected_columns = {"gene", "log2_fold_change", "p_value", "padj"}
    if not expected_columns.issubset(df.columns):
        raise ValueError(f"Differential expression table missing columns {expected_columns}")
    return df


def load_expression_long(data_dir: Optional[str | Path] = None) -> pd.DataFrame:
    base = resolve_data_dir(data_dir)
    table = base / "airway_expression_long.csv"
    if not table.exists():
        raise FileNotFoundError(f"Missing expression long table at {table}")
    df = pd.read_csv(table)
    expected_columns = {"sample_id", "gene", "log2_expression", "condition"}
    if not expected_columns.issubset(df.columns):
        raise ValueError(f"Expression long table missing columns {expected_columns}")
    return df


def load_metadata(data_dir: Optional[str | Path] = None) -> pd.DataFrame:
    base = resolve_data_dir(data_dir)
    table = base / "airway_metadata.csv"
    if not table.exists():
        raise FileNotFoundError(f"Missing metadata table at {table}")
    df = pd.read_csv(table)
    expected_columns = {"sample_id", "condition"}
    if not expected_columns.issubset(df.columns):
        raise ValueError(f"Metadata table missing columns {expected_columns}")
    return df



