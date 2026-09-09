# AQ EEG reproducibility

Analysis scripts and aggregate results accompanying **Self-reported attention to detail and visual network organization in resting-state EEG**, by Jianyi Liu, Xuan Yao, and Xiaobin Ding.

This compact release supports independent checking of the reported screening results: 20,928 tests, 320 local false-discovery-rate families, and 14 retained associations from 11 families. It also provides the original analysis and figure-generation code for methodological inspection. The leading association concerns self-reported attention to detail and mean degree of visual-labeled nodes in eyes-open beta-AEC graphs.

**No participant-level data are distributed.** The 13 CSV files contain test-level or model-level summaries, not individual observations. Their combined size is 3,215,399 bytes (about 3.07 MiB). Raw and post-ICA EEG, connectivity matrices, per-participant graph features, AQ/demographic records, design matrices, individual residuals, and the scatterplot figure are excluded.

## Verify the public results

With Python 3.10 or newer, from this repository's root:

```sh
python scripts/verify_public_results.py
python -m unittest discover -s tests -v
```

The verification uses only the Python standard library. It checks all 13 frozen CSV checksums, verifies row/sample counts and Fisher intervals, recomputes the local and broader Benjamini-Hochberg corrections, and checks arithmetic in the regression summaries. Expected output includes:

```text
tests: 20928
local_families: 320
local_hits: 14
hit_families: 11
band_pooled_survivors: 0
studywide_survivors: 0
reported_followup_units: 7
```

This is an aggregate-result audit, not re-estimation of correlations from individual data. Bootstrap, permutation, leave-one-out estimates, graph computation, and source reconstruction cannot be rerun from this public release alone.

## Contents

```text
reference_outputs/4_corr/   12 complete correlation tables and one regression table
data_manifest.json         SHA-256 hashes of the 13 aggregate CSV files
scripts/                   Public verification and summary-figure entry points
tests/                     Tests for the public verification workflow
4_corr/                    Original MATLAB statistics and residual-export scripts
2_main_process/            Original graph and post-ICA feature audit scripts
5_paper/                   Original Python figure-generation scripts
tools/                     Full-pipeline comparison helpers requiring private inputs
docs/                      Variable definitions, provenance, and access scope
```

## Figures and tables

The correlation tables support the aggregate profiles in Figures 1 and 2 and the candidate associations in Table 2. The regression table supports Table 4 and Figure 4. Figure 3 requires participant-level residuals and is not distributed; the public figure entry point never loads or renders it.

The plotting routines in `5_paper/build_figures123_v2.py` remain unchanged. Its default entry point checks private residual data and therefore cannot run on this public package alone. Use the dedicated public entry point described in `docs/REPRODUCTION_SCOPE.md` for Figures 1, 2, and 4. Historical `v2` script names refer to figure revisions; their outputs were verified against the images embedded in the v8 manuscript.

## Original MATLAB analysis

The scripts preserve the original computations. The full entry point is `run_pipeline('stats')`, but it requires the private AQ table, canonical graphs and full reference set. It is **not** the quick-start command for this public release.

The authoritative environment used MATLAB R2024b and Statistics and Machine Learning Toolbox. Graph/source workflows additionally used DISCOVER-EEG v2.1 with the true-upper-triangle density fix, BCT 2019-03-03, EEGLAB 2024.0, FieldTrip 2024-07-01, and Parallel Computing Toolbox. These third-party tools are not redistributed. Set `DISCOVER_EEG_HOME`, `BCT_HOME`, `EEGLAB_HOME`, and `FIELDTRIP_HOME` when using an authorized complete private analysis copy.

The public copy of `pipeline_config.m` removes machine-specific fallback paths. The public `tools/verify_results.m` adds explicit checks for the complete private reference set so that missing design matrices cannot produce a misleading success message. No statistical formula, screening rule, random seed, or frozen result was changed for publication.

## Interpretation of historical labels

The study was not preregistered. The machine label `family=confirm` selects the seven reported follow-up units; it does not indicate confirmatory or preregistered analysis. The eighth regression row is a historical exploratory unit and remains available for transparency, but it is not a Table 4 unit.

Likewise, `permFL` denotes the implementation-specific residual-permutation calculation described in Methods; the internal prediction procedure was conditional and restricted in scope, not a fully nested repetition of the complete discovery process. See `docs/DATA_DICTIONARY.md` before interpreting these fields.

## Access and citation

This public release contains scripts and aggregate results only. Participant-level files are not released. A procedure for requesting those files has not been established in this repository; no promise of access should be inferred.

Use `CITATION.cff` to cite this versioned software and aggregate-results companion. A GitHub release/tag and commit identify the published snapshot; no DOI has been assigned by GitHub. Original code and associated documentation are licensed under MIT; the 13 aggregate statistical CSV tables are licensed under CC BY 4.0. See [LICENSES.md](LICENSES.md) for the separate scopes, attribution and third-party exclusions. The license notices were added on 9 September 2026 and apply to the unchanged v1.0.0 materials; the release tag, scripts and data have not been altered.
