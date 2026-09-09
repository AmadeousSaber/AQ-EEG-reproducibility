# AQ EEG reproducibility

Version 1.1.0 accompanies **Self-reported attention to detail and visual network organization in resting-state EEG**, by Jianyi Liu, Xuan Yao, and Xiaobin Ding.

The repository provides analysis scripts and **aggregate results only**. Its 21 statistical CSVs comprise 13 unchanged historical tables and eight revision-specific tables (the new tables total 13,940 bytes). No participant-level AQ/demographic observations, EEG, connectivity/graph matrices, design matrices, individual residuals, deletion traces, or bootstrap draws are released.

## What changed in v1.1.0

The primary correlation screen is unchanged: 20,928 tests in 320 local false-discovery-rate families, with 14 retained associations from 11 families. The original v1.0.0 files remain available unchanged for provenance.

The seven exploratory follow-up models now use the same full-minus-reduced R² statistic on the observed outcome and every residual-permuted pseudo-outcome. Both models are refitted on each pseudo-outcome. A partial-F companion uses the identical permutations. Bootstrap intervals re-estimate standardization/PC1 within the fixed selected feature sets. Model features, observed R²/ΔR², and restricted-scope leave-one-out estimates are unchanged.

A targeted stability assessment is supplied for the eyes-open communication-by-sex interaction identified in that audit. Its deletion, robust-SE, signed-bootstrap and fixed-design null checks are **post hoc sensitivity analyses**, not independent discovery or replication. The seven-model table remains the source for the reported primary interaction p/q values; the targeted audit does not replace them with a different random-seed calculation.

## Verify without participant data

Python 3.10 or newer; no third-party packages are needed:

```sh
python scripts/verify_public_results.py
python scripts/verify_revised_results.py
python -m unittest discover -s tests -v
```

The historical checker and its 13-file `data_manifest.json` are unchanged. The revised checker validates the eight new checksums/row counts, plus-one permutation p values, separate seven-model BH corrections, interval structure, summary arithmetic and the focal stability summaries. These checks inspect published aggregates; they do not re-estimate resampling distributions from participant data.

## Contents and manuscript mapping

| Location | Contents and role |
| --- | --- |
| `reference_outputs/4_corr/` | Twelve complete correlation tables and the historical eight-row regression table; unchanged v1.0.0 provenance. |
| `reference_outputs/4_corr_v2/` | Current seven-model regression summary, 28-row paired-ΔR²/partial-F permutation table and seven-row bootstrap comparison. Source for revised Table 4 and Figure 4. |
| `reference_outputs/interaction_stability/` | Five aggregate focal-interaction diagnostics; no individual deletion trace or resampling draws. |
| `data_manifest.json`, `revised_manifest.json` | Separate historical and revision-specific SHA-256/metadata inventories. |
| `4_corr/step6_hierarchical_regression_v2.m` | Standalone revised follow-up analysis for holders of authorized private design inputs. |
| `4_corr/step7_interaction_stability.m` | Standalone targeted interaction sensitivity workflow for authorized private inputs. |
| `scripts/`, `tests/`, `docs/` | Public checking, aggregate figure entry points, tests, variable definitions and reproduction scope. |

Figures 1 and 2 use unchanged screening tables. Revised Figure 4 uses the new seven-model summary. Figure 3 depends on participant residuals and is neither distributed nor rendered by the public entry point.

## Aggregate figures

```sh
python -m pip install -r requirements-figures.txt
python scripts/render_public_figures.py --check-inputs
python scripts/render_public_figures.py
python 5_paper/build_figure4_v3.py
```

The public entry point renders current Figures 1, 2 and 4 as previews in `outputs/figures/`; `--historical` selects the historical Figure 4. The final command exports revised Figure 4 as SVG/PDF/PNG and RGB 600-dpi TIFF. It reads only the checksum-locked aggregate summary.

## MATLAB analysis and preprocessing provenance

The original MATLAB scripts remain available unchanged. `run_pipeline('stats')` invokes the historical pipeline, not the revised follow-up entry point, and requires private inputs; it is not a public quick-start command. To run the revision with an authorized complete copy:

```matlab
addpath('4_corr');
step6_hierarchical_regression_v2
step7_interaction_stability
```

Existing v2 outputs require an explicit `step6_hierarchical_regression_v2(true)` refresh. See each script for its inputs and outputs; do not publish its participant-level MAT or trace outputs.

The analysis environment used MATLAB R2024b and Statistics and Machine Learning Toolbox. Preprocessing used a **study-specific adaptation of DISCOVER-EEG**, based on an upstream distributed codebase; the historical local label “v2.1” is not an official DISCOVER-EEG release identifier. Graph construction also incorporated the verified true-upper-triangle density correction. Other external dependencies include BCT 2019-03-03, EEGLAB 2024.0 and FieldTrip 2024-07-01. Third-party tools/atlas files are not redistributed. The public `pipeline_config.m` uses portable project-relative data paths and environment variables for external software.

## Interpretation, access and license

The study was not preregistered. Historical `family=confirm` selects seven paper units but is not a confirmatory-design claim. New tables instead label them `exploratory_fixed_units`; do not apply the old `confirm` filter to the revised table. The legacy `*_permFL` column names are retained for compatibility in the revised summary, but now explicitly denote paired-ΔR² residual-permutation checks. Read [DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md) and [REPRODUCTION_SCOPE.md](docs/REPRODUCTION_SCOPE.md) for the distinct historical/current definitions.

Participant-level files are not publicly released. This repository does not establish a request/access procedure or promise later access. Full re-estimation requires authorized private inputs; there are no synthetic stand-ins. No DOI has been assigned here; cite the version and GitHub snapshot using `CITATION.cff`.

Original code/documentation are MIT licensed; all 21 aggregate CSV tables are CC BY 4.0 licensed. The original notices are retained and the new directories have explicit CC BY 4.0 notices. Third-party materials and participant-level data are outside these grants. See [LICENSES.md](LICENSES.md) and [CHANGELOG.md](CHANGELOG.md).
