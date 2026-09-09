# Change log

## 1.1.0 — 2026-09-10

- Added revised, self-contained seven-model follow-up code: paired observed/pseudo-outcome ΔR² residual permutations (10,000), partial F on identical permutations, and separate seven-main/seven-interaction BH correction.
- Re-estimated standardization/PC1 inside 5,000 fixed-feature bootstrap resamples; corrected the OLS helper's standard-error shape and replaced explicit normal-equation inversion with QR solves.
- Added three revised follow-up aggregate CSVs and five targeted interaction stability CSVs, with a separate manifest/checker/tests.
- Added targeted post hoc interaction diagnostics: delete-one perturbations, maximum-Cook deletion summary, conventional/HC3 coefficient inference, sex-stratified signed-coefficient bootstrap and iid-Gaussian fixed-design interaction-null calibration.
- Preserved the original features, observed R²/ΔR² and restricted-scope LOO results. The 13 original aggregate CSVs, historical manifest/checker and original step6 script remain unchanged.
- Added revised aggregate Figure 4 rendering and clarified current/historical figure entry points.
- Clarified study-specific DISCOVER-EEG adaptation provenance, historical field-name semantics, licenses and public/private reproduction boundaries.
- No participant-level data, residuals, MAT files, deletion traces, resampling draws or private logs were added.

The v1.0.0 tag is a historical snapshot and must not be moved. New primary p/q values come from the seven-model seed-20260910 run, not a targeted focal rerun.

## Licensing update — 2026-09-09

Added MIT notices for original code/documentation and CC BY 4.0 notices for the original 13 aggregate CSVs; no statistical data or v1.0.0 tag was changed.

## 1.0.0 — 2026-09-05

Initial scripts-and-aggregate-results release. Participant-level inputs were excluded.
