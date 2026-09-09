# Reproduction and disclosure scope

## Historical and current results

Version 1.1.0 adds eight aggregate CSV files without changing the 13 historical CSVs, their manifest, the historical verifier, or the original MATLAB analysis script. The primary 20,928-test correlation screen is unchanged.

Current Table 4/Figure 4 use `reference_outputs/4_corr_v2/results_regression_v2_summary.csv`. The old `reference_outputs/4_corr/results_regression_summary.csv` is preserved for historical reproducibility, not silently replaced. Its eight rows include a historical exploratory U0 unit; the new table contains exactly the seven paper units.

The revised permutation analysis uses seed 20260910, 10,000 permutations, observed and pseudo-outcome paired ΔR², and an identical-permutation partial-F companion. Seven main-effect and seven interaction p values are BH-adjusted separately. The revised bootstrap uses seed 20260911 and 5,000 paired resamples; standardization/PCA are re-estimated within the fixed original features. Neither analysis repeats the original outcome-guided feature discovery.

The targeted eyes-open communication interaction diagnostics use separate seeds: deletion permutations 20260913, sex-stratified signed-coefficient bootstrap 20260914, and interaction-null simulation 20260915. The primary full-sample p/q remain those in the seven-model table, not a targeted rerun with the new seed.

## Publicly reproducible checks

- Validate historical 13-file and new eight-file inventories, byte counts/rows and SHA-256 values.
- Recalculate historical local, band-pooled and study-wide BH corrections.
- Recalculate paired-permutation plus-one p values and four seven-model BH families from published exceedance counts.
- Check revised summary/companion-table links, bootstrap interval ordering, unchanged R²/LOO values and aggregate interaction diagnostic arithmetic.
- Recreate Figures 1, 2 and revised Figure 4 from aggregate values; inspect the analysis code and selected unit definitions.

The public checks validate aggregate consistency, not the correctness of unobserved participant-level input processing. Monte Carlo intervals concern simulation/permutation sampling uncertainty and do not correct for feature selection.

## Not reproducible from this release alone

Re-estimation of correlations, residual permutations, bootstrap distributions, leave-one-out predictions, EEG graph features or source reconstruction requires private participant-level inputs. Figure 3 and individual scatterplots also require private residuals.

The release excludes participant AQ/demographic records, raw/post-ICA EEG, individual connectivity/graph files, design matrices, MAT/residual tables, delete-one traces, bootstrap/simulation draws, per-person manifests, private metadata/log snapshots, the full manuscript, third-party PDFs/toolboxes, internal project memory and credential helpers. Only the five specifically named focal-interaction aggregate CSVs are public.

The authors have not publicly released participant-level inputs, and this repository does not establish an access-request mechanism. It does not claim anonymous data, controlled-access approval, consent-based sharing permissions or future access that have not been confirmed.

## Run public checks and aggregate figures

```sh
python scripts/verify_public_results.py
python scripts/verify_revised_results.py
python -m unittest discover -s tests -v
python -m pip install -r requirements-figures.txt
python scripts/render_public_figures.py --check-inputs
python scripts/render_public_figures.py
```

The public figure entry point never calls Figure 3 or reads participant residuals. Add `--historical` to use historical Figure 4; the default uses revised Table 4. PNG previews go to `outputs/figures/`. For editable/high-resolution revised Figure 4 use `python 5_paper/build_figure4_v3.py`; its source hash is locked and TIFF output is RGB/600 dpi. Exact pixels/fonts may vary by environment.

## Authorized private-input analysis

The revised `step6_hierarchical_regression_v2.m` uses `pipeline_config.m`, the atlas CSV and four private design-data MAT files; it does not depend on a private audit-work directory. It preserves the original selected features and original restricted-scope LOO computation. Its `plot_data_regression_v2.mat` contains individual residuals and must not be uploaded.

`step7_interaction_stability.m` reconstructs the existing single focal feature from authorized private design data. Deletion perturbations are diagnostics, not final-model participant exclusions or 39 additional discoveries. Its HC3 reference is an approximate t diagnostic. Its iid-Gaussian fixed-design interaction-null simulation does not assess all noise structures, discovery-wide inference or seven-interaction FDR.

The original `run_pipeline('stats')` remains historical and requires additional private AQ/graph/reference data. The public code does not supply synthetic replacements. See script headers and the data dictionary before running any full pipeline.

## Provenance and interpretive boundaries

The study was not preregistered. The original analysis freeze was 2026-08-18; revision-specific computations were audited on 2026-09-10. Preprocessing used a study-specific adaptation of DISCOVER-EEG; “v2.1” was an internal local label, not an official upstream release. External dependency licenses remain with their owners.

Revised ΔR² percentile intervals describe uncertainty within the fixed selected feature set. Because ΔR² is nonnegative for nested OLS models, exclusion of zero by such an interval is not alone a significance test. The signed focal-interaction bootstrap is also post hoc and not selection-adjusted. Sign retention is a descriptive stability proportion, not a p value.
