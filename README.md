# Gene Expression Visuals

[![CI](https://github.com/SoheiltheStar/gene-expression-visuals/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/SoheiltheStar/gene-expression-visuals/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)

Interactive RNA-seq analysis and dashboard showcasing differential expression (volcano), clustered heatmap, PCA, and a lightweight classifier. Built for portfolio/demo use: fast, reproducible, and easy to deploy.

## Dataset

- Source: GEO accession [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778) (Airway smooth muscle; Dex vs Untreated)
- Script `scripts/prepare_airway_subset.py` downloads the supplementary FPKM matrix and generates tidy processed tables under `data/processed/`.
- Processed artifacts are included for immediate exploration; re-run the script to regenerate as needed.

## Quickstart

```bash
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .[dev]
# optional: refresh processed data
python scripts/prepare_airway_subset.py
streamlit run app/streamlit_app.py
```

Optional: `pytest` to run smoke tests.

## Visual Story Highlights

- Interactive volcano plot with significance and effect filters
- Clustered heatmap of the top differentially expressed genes
- PCA scatter with condition-aware shading and sample hovercards
- Logistic regression classifier with feature importances and ROC curve

## Project Layout

```
gene-expression-visuals/
├── app/                     # Streamlit UI entrypoint
├── data/
│   ├── processed/           # Tidy subsets for visuals (generated)
│   └── raw/                 # GEO downloads (gitignored)
├── docs/                    # Narrative notes and visuals exports
├── notebooks/               # Optional exploratory notebooks (templates)
├── scripts/                 # CLI utilities (data prep)
├── src/gene_expression_visuals/
│   ├── __init__.py
│   ├── analysis.py
│   ├── data_loader.py
│   ├── ml.py
│   └── visualization.py
├── tests/                   # Pytest-based smoke checks
├── pyproject.toml
├── README.md
└── .env.example
```

## Deployment

- Streamlit Cloud: push this repo and set the app entrypoint to `app/streamlit_app.py`. No secrets required.
- Local container (optional): create a minimal Dockerfile if needed; Streamlit runs on port 8501.

## Portfolio Positioning

- End-to-end ownership: data ingest → stats → ML → interactive storytelling
- Reproducible and maintainable: scripted prep, tests, packaging, and CI
- Easy to extend: survival curves, enrichment, or additional biomarker panels

Repo: `https://github.com/SoheiltheStar/gene-expression-visuals`

