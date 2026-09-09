"""Verify v1.1.0 aggregate tables; never load participant-level data.

This standard-library audit verifies hashes and arithmetic, not resampling
from private observations. Output contains summary counts only.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

from verify_public_results import bh_adjust, require


ROOT = Path(__file__).resolve().parents[1]
PREFIX = "reference_outputs/4_corr_v2/"
STABILITY = "reference_outputs/interaction_stability/"
EXPECTED_ROWS = {
    PREFIX + "results_regression_v2_summary.csv": 7,
    PREFIX + "results_regression_v2_permutations.csv": 28,
    PREFIX + "results_regression_v2_bootstrap.csv": 7,
    STABILITY + "observed_interaction.csv": 1,
    STABILITY + "delete_one_summary.csv": 1,
    STABILITY + "maximum_cook_deletion.csv": 1,
    STABILITY + "bootstrap_summary.csv": 1,
    STABILITY + "interaction_null_summary.csv": 3,
}
EXPECTED_UNITS = {
    ("open", "detail"), ("open", "switching"), ("open", "communication"),
    ("closed", "detail"), ("closed", "switching"),
    ("closed", "communication"), ("closed", "imagination"),
}
METHODS = {"paired_delta_R2", "partial_F"}


def number(row: dict, key: str) -> float:
    value = float(row[key])
    require(math.isfinite(value), f"Nonfinite field: {key}")
    return value


def equal(left: float, right: float, label: str, tolerance: float = 1e-11) -> None:
    require(math.isfinite(left) and math.isfinite(right), f"Nonfinite {label}")
    require(abs(left - right) <= tolerance, f"Arithmetic mismatch: {label}")


def integer(row: dict, key: str, minimum: int = 0) -> int:
    value = number(row, key)
    require(value.is_integer() and value >= minimum, f"Invalid count: {key}")
    return int(value)


def plus_one(count: int, permutations: int) -> float:
    require(isinstance(count, int) and isinstance(permutations, int), "Counts must be integers")
    require(permutations > 0 and 0 <= count <= permutations, "Invalid permutation counts")
    return (count + 1) / (permutations + 1)


def wilson(count: int, trials: int) -> tuple[float, float]:
    require(trials > 0 and 0 <= count <= trials, "Invalid binomial counts")
    z = 1.959963984540054
    p = count / trials
    denominator = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / denominator
    half = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def indexed(rows: list[dict], columns: tuple[str, ...]) -> dict:
    result = {tuple(row[column] for column in columns): row for row in rows}
    require(len(result) == len(rows), "Duplicate aggregate row keys")
    return result


def read_tables(root: Path = ROOT) -> dict[str, list[dict]]:
    manifest = json.loads((root / "revised_manifest.json").read_text(encoding="utf-8"))
    require(manifest["version"] == "1.1.0", "Unexpected revised version")
    require(manifest["date"] == "2026-09-10", "Unexpected revision date")
    require(manifest["participant_level_data_included"] is False, "Unexpected access scope")
    require(manifest["license"] == "CC-BY-4.0", "Unexpected aggregate license")
    require(manifest["seeds"] == {
        "seven_unit_permutations": 20260910, "seven_unit_bootstrap": 20260911,
        "focal_deletion_permutations": 20260913, "focal_stratified_bootstrap": 20260914,
        "focal_interaction_null": 20260915,
    }, "Unexpected random-seed metadata")
    require(manifest["resampling"] == {
        "seven_unit_permutations": 10000, "seven_unit_bootstrap": 5000,
        "focal_deletion_permutations": 10000, "focal_stratified_bootstrap": 5000,
        "focal_null_simulations": 1000, "focal_null_inner_permutations": 999,
    }, "Unexpected resampling metadata")
    entries = manifest["files"]
    mapped = {entry["path"]: entry for entry in entries}
    require(len(mapped) == len(entries) == 8 and set(mapped) == set(EXPECTED_ROWS),
            "Unexpected revised manifest paths")
    # Exact relative allowlist above prevents manifest path traversal.
    for folder in (PREFIX, STABILITY):
        actual = {path.relative_to(root).as_posix() for path in (root / folder).rglob("*.csv")}
        expected = {path for path in EXPECTED_ROWS if path.startswith(folder)}
        require(actual == expected, "Unexpected CSV inventory in revised directories")
        require(not any((root / folder).rglob("*.mat")), "Participant MAT file in aggregate directory")
    tables = {}
    for path, row_count in EXPECTED_ROWS.items():
        entry = mapped[path]
        content = (root / path).read_bytes()
        require(len(content) == entry["bytes"], f"Byte-count mismatch: {path}")
        require(hashlib.sha256(content).hexdigest() == entry["sha256"], f"Checksum mismatch: {path}")
        with (root / path).open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        require(len(rows) == entry["rows"] == row_count, f"Row-count mismatch: {path}")
        tables[path] = rows
    return tables


def check_tables(tables: dict[str, list[dict]]) -> dict:
    require(set(tables) == set(EXPECTED_ROWS), "Unexpected table keys")
    for path, expected in EXPECTED_ROWS.items():
        require(len(tables[path]) == expected, "Unexpected aggregate row count")
    summary = indexed(tables[PREFIX + "results_regression_v2_summary.csv"], ("cond", "aqvar"))
    require(set(summary) == EXPECTED_UNITS, "Unexpected revised model set")
    permutations = indexed(tables[PREFIX + "results_regression_v2_permutations.csv"],
                           ("id", "effect", "method"))
    expected_permutation_keys = {
        (condition + "_" + domain, effect, method)
        for condition, domain in EXPECTED_UNITS
        for effect in ("main", "interaction") for method in METHODS
    }
    require(set(permutations) == expected_permutation_keys, "Unexpected permutation groups")
    for key, row in permutations.items():
        count, trials = integer(row, "exceedances"), integer(row, "B", 1)
        require(trials == 10000, "Unexpected seven-unit permutation count")
        equal(number(row, "p_plus_one"), plus_one(count, trials), "plus-one p")
        lo, hi = wilson(count, trials)
        equal(number(row, "tail_probability_wilson95_lo"), lo, "permutation Wilson lower")
        equal(number(row, "tail_probability_wilson95_hi"), hi, "permutation Wilson upper")
        require(number(row, "observed_statistic") >= 0, "Negative observed nested statistic")
        require(row["id"] == row["cond"] + "_" + row["aqvar"], "Permutation identity mismatch")
        require(integer(row, "n") == (39 if row["cond"] == "open" else 41), "Permutation n mismatch")
    for effect in ("main", "interaction"):
        for method in METHODS:
            group = [row for key, row in permutations.items() if key[1:] == (effect, method)]
            require(len(group) == 7, "BH family is not seven models")
            for row, adjusted in zip(group, bh_adjust([number(row, "p_plus_one") for row in group]), strict=True):
                equal(number(row, "q_BH_seven"), adjusted, "seven-model BH")

    bootstrap = indexed(tables[PREFIX + "results_regression_v2_bootstrap.csv"], ("id",))
    require(set(bootstrap) == {(condition + "_" + domain,) for condition, domain in EXPECTED_UNITS},
            "Unexpected bootstrap model set")
    for (condition, domain), row in summary.items():
        identifier = condition + "_" + domain
        require(row["family"] == "exploratory_fixed_units", "Historical family filter applied to revised data")
        require(integer(row, "N") == (39 if condition == "open" else 41), "Revised n mismatch")
        require(integer(row, "B_permutations") == 10000 and integer(row, "B_bootstrap") == 5000,
                "Summary resampling counts differ")
        equal(number(row, "dR2_main"), number(row, "R2_main") - number(row, "R2_base"), "main delta R2")
        equal(number(row, "dR2_inter"), number(row, "R2_inter") - number(row, "R2_main"), "interaction delta R2")
        equal(number(row, "dq2"), number(row, "q2_unit") - number(row, "q2_base"), "LOO difference")
        for effect, method, p_column, q_column in (
            ("main", "paired_delta_R2", "dR2_main_p_permFL", "q_unit"),
            ("interaction", "paired_delta_R2", "inter_p_permFL", "q_inter"),
            ("main", "partial_F", "dR2_main_p_permF", "q_unit_F"),
            ("interaction", "partial_F", "inter_p_permF", "q_inter_F"),
        ):
            check = permutations[(identifier, effect, method)]
            equal(number(row, p_column), number(check, "p_plus_one"), "summary p linkage")
            equal(number(row, q_column), number(check, "q_BH_seven"), "summary q linkage")
            if method == "paired_delta_R2":
                equal(number(check, "observed_statistic"), number(row, "dR2_main" if effect == "main" else "dR2_inter"),
                      "observed delta linkage")
        boot = bootstrap[(identifier,)]
        require(integer(boot, "Bboot") == 5000 and integer(boot, "paired_valid_draws") == 5000,
                "Expected 5,000 valid paired bootstraps in this snapshot")
        equal(number(row, "bootstrap_valid_draws"), number(boot, "paired_valid_draws"), "bootstrap count linkage")
        for prefix in ("fixed_component", "refit_component"):
            require(0 <= number(boot, prefix + "_ci_lo") <= number(boot, prefix + "_ci_hi") <= 1,
                    "Invalid bootstrap interval")
        equal(number(row, "dR2_boot_lo"), number(boot, "refit_component_ci_lo"), "bootstrap lower linkage")
        equal(number(row, "dR2_boot_hi"), number(boot, "refit_component_ci_hi"), "bootstrap upper linkage")

    # This is a version-specific snapshot check, not an inferential threshold.
    # It prevents substitution of the separate seed-20260913 focal rerun.
    focal = summary[("open", "communication")]
    equal(number(focal, "inter_p_permFL"), plus_one(41, 10000), "primary focal paired p")
    equal(number(focal, "inter_p_permF"), plus_one(40, 10000), "primary focal F p")
    observed = tables[STABILITY + "observed_interaction.csv"][0]
    require(integer(observed, "n") == 39 and integer(observed, "male_n") + integer(observed, "female_n") == 39,
            "Focal sample counts differ")
    equal(number(observed, "delta_R2_interaction"), number(focal, "dR2_inter"), "focal delta linkage")
    equal(number(observed, "R2_full") - number(observed, "R2_reduced"), number(observed, "delta_R2_interaction"),
          "focal nested-model arithmetic")
    beta = number(observed, "beta_interaction_AQ_points_per_brain_SD")
    require(-1.7 < beta < -1.5 and .14 < number(observed, "delta_R2_interaction") < .15,
            "Unexpected focal coefficient/variance snapshot")
    for prefix in ("conventional", "HC3"):
        equal(number(observed, prefix + "_t"), beta / number(observed, prefix + "_SE"), "coefficient t statistic")
        require(number(observed, prefix + "_CI_lo") < beta < number(observed, prefix + "_CI_hi") < 0,
                "Unexpected signed interaction interval")
    equal(number(observed, "partial_F"), number(observed, "conventional_t") ** 2, "one-df F/t identity")

    deletion = tables[STABILITY + "delete_one_summary.csv"][0]
    require(integer(deletion, "deletions") == integer(deletion, "same_direction_count") == 39,
            "Unexpected deletion/sign-retention summary")
    for prefix in ("beta", "delta_R2", "paired_delta_p", "partial_F_p"):
        require(number(deletion, prefix + "_min") <= number(deletion, prefix + "_median") <= number(deletion, prefix + "_max"),
                "Invalid deletion range")
    for method in ("paired_delta", "partial_F"):
        require(0 <= number(deletion, method + "_p_min") <= number(deletion, method + "_p_max") <= 1,
                "Invalid deletion p range")
        require(integer(deletion, method + "_p_le_05_count") <= 39, "Invalid threshold count")
    maximum = tables[STABILITY + "maximum_cook_deletion.csv"][0]
    require(integer(maximum, "retained_n") == 38 and integer(maximum, "interaction_sign_retained") == 1,
            "Invalid maximum-Cook perturbation")
    equal(number(maximum, "original_Cook_D"), number(observed, "max_Cook_D"), "Cook summary linkage")
    for method in ("paired_delta", "partial_F"):
        count = integer(maximum, method + "_exceedances")
        equal(number(maximum, method + "_p"), plus_one(count, 10000), "deletion plus-one p")
        lo, hi = wilson(count, 10000)
        equal(number(maximum, method + "_tail_Wilson95_lo"), lo, "deletion Wilson lower")
        equal(number(maximum, method + "_tail_Wilson95_hi"), hi, "deletion Wilson upper")
    signed = tables[STABILITY + "bootstrap_summary.csv"][0]
    require(integer(signed, "bootstrap_draws") == integer(signed, "valid_draws") == 5000,
            "Invalid signed-bootstrap count")
    equal(number(signed, "observed_beta_brain_standardized"), beta, "bootstrap observed beta")
    require(number(signed, "percentile95_beta_lo") < beta < number(signed, "percentile95_beta_hi") < 0,
            "Invalid signed-bootstrap interval")
    equal(number(signed, "fraction_beta_positive") + number(signed, "fraction_beta_negative"), 1,
          "bootstrap sign fractions")
    equal(number(signed, "descriptive_sign_retention"), number(signed, "fraction_beta_negative"),
          "descriptive sign retention")

    null_rows = indexed(tables[STABILITY + "interaction_null_summary.csv"], ("method",))
    require(set(null_rows) == {(method,) for method in METHODS | {"original_constant_reduced_R2"}},
            "Unexpected null-calibration methods")
    for row in null_rows.values():
        count, simulations = integer(row, "rejections"), integer(row, "simulations", 1)
        require(simulations == 1000 and integer(row, "inner_permutations") == 999, "Unexpected null design counts")
        equal(number(row, "rejection_rate"), count / simulations, "null rejection rate")
        lo, hi = wilson(count, simulations)
        equal(number(row, "Wilson95_lo"), lo, "null Wilson lower")
        equal(number(row, "Wilson95_hi"), hi, "null Wilson upper")
        require(number(row, "generating_residual_SD") > 0, "Invalid generating SD")
    return {"status": "PASS", "revised_csv_files": 8, "revised_followup_units": 7,
            "permutation_effect_statistic_rows": 28, "bootstrap_units": 7,
            "focal_stability_files": 5, "participant_level_inputs_loaded": 0}


def audit(root: Path = ROOT) -> dict:
    tables = read_tables(root)
    report = check_tables(tables)
    historical_path = root / "reference_outputs" / "4_corr" / "results_regression_summary.csv"
    historical_manifest = json.loads((root / "data_manifest.json").read_text(encoding="utf-8"))
    require(hashlib.sha256(historical_path.read_bytes()).hexdigest() ==
            historical_manifest["results_regression_summary.csv"], "Historical regression checksum mismatch")
    with historical_path.open(encoding="utf-8-sig", newline="") as handle:
        historical = indexed([row for row in csv.DictReader(handle) if row["family"] == "confirm"],
                             ("cond", "aqvar"))
    require(set(historical) == EXPECTED_UNITS, "Unexpected historical model keys")
    current = indexed(tables[PREFIX + "results_regression_v2_summary.csv"], ("cond", "aqvar"))
    for key, row in current.items():
        for column in ("N", "R2_base", "R2_main", "R2_inter", "dR2_main", "dR2_inter", "std_beta",
                       "q2_base", "q2_unit", "dq2"):
            equal(number(row, column), number(historical[key], column), "historical linkage " + column)
    report["historical_effect_and_loo_linkage"] = "PASS"
    return report


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
