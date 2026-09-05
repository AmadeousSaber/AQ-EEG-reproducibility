# Reproduction and disclosure scope

## Available from this release

- Check the exact inventory, SHA-256 values and structure of the 13 frozen aggregate CSV files.
- Recalculate local, band-pooled and study-wide BH adjustments from all 20,928 published p values.
- Check confidence-interval arithmetic, the visual-network focal estimate and conditional-prediction difference arithmetic.
- Recreate aggregate Figures 1, 2 and 4 from the published summary tables.
- Inspect the MATLAB scripts used for feature extraction, correlations, local correction, residual exports, and selection-conditioned follow-up models.

## Not available from this release alone

Re-estimating correlations, bootstrap/permutation distributions, leave-one-out models, graph measures, source reconstruction or individual scatterplots requires participant-level inputs that the authors have chosen not to release publicly. The repository does not contain synthetic stand-ins and does not claim end-to-end public reproduction.

## Summary-figure entry point

Install the tested plotting dependencies in your preferred Python environment, then run:

```sh
python -m pip install -r requirements-figures.txt
python scripts/render_public_figures.py --check-inputs
python scripts/render_public_figures.py
```

The input-check command validates public inputs without drawing. The final command invokes the existing Figure 1, 2 and 4 plotting functions and writes PNG previews to `outputs/figures/`. It never calls Figure 3 or reads `plot_data_all.csv`. Plotted numerical inputs match the manuscript; exact pixels and font metrics can vary by software/font environment. Python 3.12.14 with the versions in `requirements-figures.txt` passed the public input check.

## Intentional exclusions

The release whitelist excludes all participant-level AQ and demographic tables, raw/post-ICA EEG, individual connectivity and graph files, design matrices, residual tables/MAT files, individual manifests, parameter/log snapshots, the full manuscript and its Figure 3 scatterplot, third-party PDFs/toolboxes, internal project memory, and GitHub credential/diagnostic helpers.

The four uncorrected-significance subset tables are omitted to avoid duplication; they can be recovered by filtering the complete nodal tables at `p < .05` and are not the manuscript's corrected evidence set.

## Changes made solely for distribution

- Removed machine-specific external-toolbox fallback paths from `pipeline_config.m`; environment-variable configuration is retained.
- Added fail-closed checks to `tools/verify_results.m` for the private four-design-matrix/eighteen-CSV reference set, avoiding incomplete-verification success.
- Added public aggregate auditing/tests and an aggregate-only figure entry point.
- Updated figure-script usage docstrings to point to the tested public entry point; plotting routines are unchanged.
- No original statistical formula, graph threshold, selection rule, random seed, correlation value or regression estimate was altered.

The analysis freeze date is 2026-08-18. The study was not preregistered. Code labels such as `confirm`, `permFL`, and `nestedLOO` must be read using the manuscript-aligned definitions in `DATA_DICTIONARY.md`.
