# Gene Expression Insights — Case Study

## Challenge

Quantify how airway smooth muscle cells respond to dexamethasone therapy compared with untreated controls, then package the insight in a decision-friendly narrative.

## Approach

- Data ingest: downloaded the GEO GSE52778 FPKM matrix (Dex vs Untreated) via `scripts/prepare_airway_subset.py` for reproducible preprocessing.
- Signal enrichment: log transformed the expression matrix and ran Welch t-tests with Benjamini–Hochberg correction to prioritise significant genes.
- Storytelling: delivered interactive volcano, heatmap, PCA, and classifier diagnostics through a Streamlit experience.

## Highlights

- >1,000 genes reach FDR < 0.05; FKBP5, PER1, and GILZ headline the response signature.
- PCA separates treated vs control replicates cleanly, confirming strong transcriptional shifts.
- Logistic regression on the top 25 genes scores ROC AUC > 0.85, indicating deployable biomarker potential.

## Portfolio Usage

- Embed dashboard screenshots or GIFs to illustrate interactivity.
- Reference reproducibility (scripted download, tests, pyproject packaging) in proposal copy.
- Tie the narrative to pharma/biotech needs: steroid response assessment, biomarker discovery, mechanism validation.

