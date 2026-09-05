# Aggregate result dictionary

## Observation unit and provenance

Each correlation-table row is one age- and sex-adjusted association across 39 eyes-open or 41 eyes-closed participants. Each regression-table row is a model-level summary. No row is a participant record. CSV files are byte-identical to the frozen reference outputs used for the manuscript.

The 12 correlation tables contain 384 global, 1,344 network-mean and 19,200 nodal tests. Files encode condition (`open`, `closed`) and estimator (`aec`, `dwpli`) in their names. AEC denotes orthogonalized amplitude-envelope correlation; dwPLI denotes debiased weighted phase-lag index.

## Correlation fields

| Field | Meaning |
| --- | --- |
| `metric` | Global metric key, including estimator, band and graph quantity. |
| `quantity` | Nodal graph quantity: `auc_degree` or `auc_cc`. |
| `band` | theta (4–7.9 Hz), alpha (8–12.9 Hz), beta (13–30 Hz), gamma (30.1–80 Hz). |
| `network` | Schaefer-7 grouping: Vis, SomMot, DorsAttn, SalVentAttn, Limbic, Cont, Default. |
| `node`, `roi` | Atlas-label index and name; these identify cortical parcels, not people. |
| `aqvar` | `total`, `social`, `switching`, `detail`, `communication`, or `imagination`; self-report content scores. |
| `r` | Pearson partial correlation after residualizing age and sex. |
| `ci_lo`, `ci_hi` | Approximate 95% Fisher interval, using two covariates. |
| `p` | Two-sided uncorrected partial-correlation p value. |
| `n` | Number of participants contributing to the test, not a participant identifier. |
| `q4` | Historical name for band-separated local BH-adjusted p value. It is not correction across all four bands. |
| `level`, `cond_conn` | Spatial scale and, where present, condition-estimator label. |

Graph measures were summarized across true graph densities of 15%, 20%, 25%, and 30%. Global quantities are clustering (`auc_gcc`), efficiency (`auc_geff`), characteristic path length (`auc_cpl`), and small-worldness (`auc_sw`). Degree is measured on binary fixed-density graphs; it is not absolute AEC strength. Network means average whole-graph nodal measures over member atlas labels, not within-network subgraphs.

## Multiple-comparison families

| Scale | One band-separated family | Tests in a family |
| --- | --- | --- |
| Global | condition × estimator × metric × band | 6 AQ scores |
| Network mean | condition × estimator × nodal quantity × band | 7 networks × 6 AQ scores = 42 |
| Node | condition × estimator × nodal quantity × network × band | parcels in that network × 6 AQ scores |

The 320 local families are the primary reporting units. Pooling four bands within each otherwise corresponding family gives 80 broader families; study-wide correction pools all 20,928 tests. Local q values do not establish study-wide significance.

## Regression fields

`results_regression_summary.csv` has eight rows. Select `family=confirm` for the seven units in Table 4/Figure 4; the label is historical, not a preregistration claim. Other rows are not additional confirmatory evidence.

- `cond`, `unit`, `family`, `aqvar`, `N`, `brain_component`: model descriptor, target and sample size.
- `R2_base`, `R2_main`, `dR2_main`, `std_beta`: selection-conditioned in-sample fit and standardized coefficient.
- `dR2_main_p_asy`, `dR2_main_p_permFL`, `dR2_boot_lo`, `dR2_boot_hi`: asymptotic and implementation-specific residual-permutation checks and bootstrap interval. Screening was outcome-guided; these are not independent confirmation.
- `R2_inter`, `dR2_inter`, `inter_p_asy`, `inter_p_permFL`: added sex-by-brain interaction fit and checks.
- `slope_male`, `slope_male_p`, `slope_female`, `slope_female_p`: descriptive group-slope estimates, not individual records or independent subgroup findings.
- `maxCookD`, `dR2_excl_maxD`, `dR2_excl_maxD_p`: aggregate influence/sensitivity summaries. No excluded participant identifier is included.
- `q2_base`, `q2_unit`, `dq2`: conditional restricted-scope leave-one-out prediction relative to the sample-mean benchmark, for baseline and brain-augmented models, and their difference. Positive `dq2` does not by itself imply positive `q2_unit`.
- `q_unit`, `q_inter`: BH-adjusted permutation checks for the seven main effects and interactions, respectively.

`NaN` denotes unavailable/not-estimable entries where present. The public checker requires finite values for every correlation estimate and its p/q/CI/n fields. The full participant-level resampling calculations are available as source code but cannot be rerun from these aggregate values.
