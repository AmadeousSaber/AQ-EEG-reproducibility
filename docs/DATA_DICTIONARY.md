# Aggregate result dictionary

## Observation unit and provenance

Each correlation-table row is one age- and sex-adjusted association across 39 eyes-open or 41 eyes-closed usable records. Each regression-table row is a model-level summary. No row is a participant record. The original 13 CSV files are byte-identical to their frozen historical references; eight revision-specific tables are separately inventoried in `revised_manifest.json`.

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

## Historical regression fields (`reference_outputs/4_corr/`)

`results_regression_summary.csv` has eight rows. Select `family=confirm` for the seven historical paper units; the label is historical, not a preregistration claim. Other rows are not additional confirmatory evidence. This historical table is not the revised Table 4/Figure 4 source.

- `cond`, `unit`, `family`, `aqvar`, `N`, `brain_component`: model descriptor, target and sample size.
- `R2_base`, `R2_main`, `dR2_main`, `std_beta`: selection-conditioned in-sample fit and standardized coefficient.
- `dR2_main_p_asy`, `dR2_main_p_permFL`, `dR2_boot_lo`, `dR2_boot_hi`: asymptotic and implementation-specific residual-permutation checks and bootstrap interval. Screening was outcome-guided; these are not independent confirmation.
- `R2_inter`, `dR2_inter`, `inter_p_asy`, `inter_p_permFL`: added sex-by-brain interaction fit and checks.
- `slope_male`, `slope_male_p`, `slope_female`, `slope_female_p`: descriptive group-slope estimates, not individual records or independent subgroup findings.
- `maxCookD`, `dR2_excl_maxD`, `dR2_excl_maxD_p`: aggregate influence/sensitivity summaries. No excluded participant identifier is included.
- `q2_base`, `q2_unit`, `dq2`: conditional restricted-scope leave-one-out prediction relative to the sample-mean benchmark, for baseline and brain-augmented models, and their difference. Positive `dq2` does not by itself imply positive `q2_unit`.
- `q_unit`, `q_inter`: BH-adjusted permutation checks for the seven main effects and interactions, respectively.

`NaN` denotes unavailable/not-estimable entries where present. The public checker requires finite values for every correlation estimate and its p/q/CI/n fields. The full participant-level resampling calculations are available as source code but cannot be rerun from these aggregate values.

## Revised follow-up fields (`reference_outputs/4_corr_v2/`)

`results_regression_v2_summary.csv` contains exactly seven selected exploratory units: EO detail, switching and communication, then EC detail, switching, communication and imagination. `family=exploratory_fixed_units` replaces the old machine label; filtering on `confirm` would incorrectly remove all revised rows.

Observed R²/ΔR², standardized coefficients and original conditional LOO values are unchanged. The following compatibility-named fields have revised, explicit meanings:

| Fields | Revised definition |
| --- | --- |
| `dR2_main_p_permFL`, `inter_p_permFL` | Plus-one residual-permutation p values using paired full-minus-reduced R² on every pseudo-outcome. Both models are refitted on that same outcome. |
| `q_unit`, `q_inter` | Separate BH corrections across the seven main effects and seven interactions using the preceding paired-ΔR² p values. |
| `dR2_main_p_permF`, `inter_p_permF`, `q_unit_F`, `q_inter_F` | Companion partial-F permutation p/q values from the identical pseudo-outcomes, not an additional selected method. |
| `dR2_boot_lo`, `dR2_boot_hi` | Percentile 95% interval after standardization/PC1 are re-estimated in each bootstrap, retaining the same original features. |
| `B_permutations`, `B_bootstrap`, `bootstrap_valid_draws` | 10,000 permutations; 5,000 attempted bootstraps; number of valid paired bootstrap draws. |

These remain selection-conditioned follow-ups. Because nested-model ΔR² is nonnegative, its percentile interval excluding zero is not by itself a valid zero-effect test. It does not account for feature discovery. `q2_base`, `q2_unit` and `dq2` retain the full-sample-mean denominator used originally; positive improvement does not imply positive absolute q² or independent prediction.

`results_regression_v2_permutations.csv` has 28 rows: seven units × two effects × two statistics. `observed_statistic` is ΔR² or partial F according to `method`. `exceedances=k`, `B=10000`, and `p_plus_one=(k+1)/(B+1)`. `tail_probability_wilson95_lo/hi` describe Monte Carlo uncertainty in the underlying exceedance probability k/B, not selection-adjusted uncertainty. `q_BH_seven` adjusts only the corresponding effect/statistic's seven p values. Seed: 20260910.

`results_regression_v2_bootstrap.csv` has seven rows. It compares intervals for an observed fixed brain component with intervals re-estimating standardization/PCA inside the same resamples. `paired_valid_draws` requires both fits to be valid; degeneracy/rank/zero-scale/tied-PC counters are explicit. `mean_paired_delta_difference` and `max_abs_paired_delta_difference` compare re-estimated minus fixed component ΔR². Negative-draw counts use a 1e−12 numerical tolerance. Seed: 20260911. The revised summary reports the re-estimated interval.

## Focal interaction stability (`reference_outputs/interaction_stability/`)

These five files summarize a targeted post hoc check of eyes-open AQ communication and the existing node-42 dwPLI-alpha clustering feature. Sex is coded 0=male, 1=female. The original brain metric is standardized without outcome-directed sign reversal. The signed coefficient is the **female-minus-male brain slope difference in AQ points per brain SD**; dividing it by AQ SD changes its scale, not the binary sex coding. A negative interaction does not imply a causal sex mechanism or independent subgroup discovery.

| File | Contents |
| --- | --- |
| `observed_interaction.csv` | One full-sample aggregate row: n=39 (19 male/20 female), nested R² and ΔR², partial F, signed interaction coefficient, conventional and HC3 SE/t/p/95% intervals, df and leverage/Cook summaries. HC3 uses an approximate t reference. |
| `delete_one_summary.csv` | One aggregate row across 39 delete-one perturbations: coefficient/ΔR²/p ranges and medians, sign retention and unadjusted p-threshold counts. These are diagnostics, not 39 independent discoveries or exclusions. |
| `maximum_cook_deletion.csv` | One anonymous maximum-Cook perturbation summary, retained n=38, coefficients, ΔR²/F and paired-ΔR²/partial-F permutation counts/p/MC intervals. No identity or row index is given. B=10000; seed 20260913 covers the focal full-sample then all deletion runs. |
| `bootstrap_summary.csv` | One row from 5,000 sex-stratified resamples retaining 19/20 group sizes: signed-coefficient intervals, valid/degenerate counts and descriptive sign retention. Seed 20260914; sign retention is not a p value. |
| `interaction_null_summary.csv` | Three statistic rows from 1,000 iid-Gaussian synthetic interaction-null datasets with the observed age/sex/brain main effects retained; 999 inner permutations, rejection counts/rates at .05, Wilson intervals and generating residual SD. Seed 20260915. |

The interaction-null exercise assesses this selected fixed design under the specified Gaussian noise only, not discovery-wide error control, seven-interaction FDR or all error distributions. Primary manuscript p/q values come from the seven-unit seed-20260910 table. The targeted full-sample rerun and private deletion/bootstrap/simulation traces are deliberately not substituted for those primary values or distributed.
