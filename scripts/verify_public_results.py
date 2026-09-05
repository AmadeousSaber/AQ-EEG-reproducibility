"""Audit the manuscript's aggregate screening results without participant data.

Uses only the Python standard library. This checks reported estimates and
recomputes multiplicity adjustments; it does not re-estimate correlations.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reference_outputs" / "4_corr"
CONDITIONS = ("open", "closed")
ESTIMATORS = ("aec", "dwpli")
LEVELS = ("global", "network", "node")
NETWORKS = {"Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def bh_adjust(values: list[float]) -> list[float]:
    require(all(math.isfinite(p) and 0 <= p <= 1 for p in values), "Invalid p value")
    count = len(values)
    order = sorted(range(count), key=values.__getitem__)
    adjusted = [0.0] * count
    ceiling = 1.0
    for rank in range(count, 0, -1):
        index = order[rank - 1]
        ceiling = min(ceiling, values[index] * count / rank)
        adjusted[index] = ceiling
    return adjusted


def verify_data_hashes() -> None:
    manifest = json.loads((ROOT / "data_manifest.json").read_text(encoding="utf-8"))
    expected_names = {f"results_{level}_{cond}_{conn}.csv"
                      for level in LEVELS for cond in CONDITIONS for conn in ESTIMATORS}
    expected_names.add("results_regression_summary.csv")
    require(set(manifest) == expected_names, "Unexpected aggregate manifest entries")
    actual = {path.name for path in DATA.glob("*.csv")}
    require(actual == expected_names, "Unexpected or missing public CSV files")
    for name, expected in manifest.items():
        require(sha256((DATA / name).read_bytes()).hexdigest() == expected,
                f"Changed frozen aggregate file: {name}")


def load_correlations() -> list[dict]:
    rows = []
    for level in LEVELS:
        for condition in CONDITIONS:
            for estimator in ESTIMATORS:
                path = DATA / f"results_{level}_{condition}_{estimator}.csv"
                with path.open(encoding="utf-8-sig", newline="") as stream:
                    table = list(csv.DictReader(stream))
                require(len(table) == {"global": 96, "network": 336, "node": 4800}[level],
                        f"Unexpected row count: {path.name}")
                for row in table:
                    require(row["level"] == level, "Scale label mismatch")
                    row.update(condition=condition, estimator=estimator)
                    for key in ("r", "ci_lo", "ci_hi", "p", "q4", "n"):
                        row[key] = float(row[key])
                        require(math.isfinite(row[key]), f"Non-finite {key}: {path.name}")
                    require(row["n"] == (39 if condition == "open" else 41), "Sample size mismatch")
                    require(-1 <= row["ci_lo"] <= row["r"] <= row["ci_hi"] <= 1,
                            "Invalid correlation interval")
                    if level == "node":
                        network = row["roi"].split("_")[2]
                        require(network in NETWORKS, "Unknown atlas-network label")
                        row["network"] = network
                    rows.append(row)
    return rows


def family_key(row: dict, pool_bands: bool = False) -> tuple:
    key = (row["condition"], row["estimator"], row["level"])
    if row["level"] == "global":
        prefix = f"{row['estimator']}_{row['band']}_"
        require(row["metric"].startswith(prefix), "Unexpected global metric name")
        key += (row["metric"][len(prefix):],)
    else:
        key += (row["quantity"],)
        if row["level"] == "node":
            key += (row["network"],)
    return key if pool_bands else key + (row["band"],)


def grouped_adjustments(rows: list[dict], pool_bands: bool = False) -> tuple[list[float], int]:
    groups = defaultdict(list)
    for index, row in enumerate(rows):
        groups[family_key(row, pool_bands)].append(index)
    result = [0.0] * len(rows)
    for indices in groups.values():
        adjusted = bh_adjust([rows[index]["p"] for index in indices])
        for index, value in zip(indices, adjusted, strict=True):
            result[index] = value
    return result, len(groups)


def audit() -> dict:
    verify_data_hashes()
    rows = load_correlations()
    local, families = grouped_adjustments(rows)
    error = max(abs(value - row["q4"]) for row, value in zip(rows, local, strict=True))
    require(error < 1e-12, "Local BH values differ from the frozen results")
    hits = [row for row, value in zip(rows, local, strict=True) if value < .05]
    hit_families = len({family_key(row) for row in hits})
    pooled, pooled_families = grouped_adjustments(rows, pool_bands=True)
    studywide = bh_adjust([row["p"] for row in rows])
    require((len(rows), families, len(hits), hit_families) == (20928, 320, 14, 11),
            "Primary manuscript screening totals differ")
    require(pooled_families == 80 and not any(q < .05 for q in pooled + studywide),
            "Broader correction result differs")
    for row in rows:
        z = math.atanh(max(-.999999, min(.999999, row["r"])))
        delta = 1.96 / math.sqrt(row["n"] - 2 - 3)
        require(abs(math.tanh(z - delta) - row["ci_lo"]) < 1e-12
                and abs(math.tanh(z + delta) - row["ci_hi"]) < 1e-12,
                "Fisher interval differs from the documented two-covariate formula")
    focal = [row for row in hits if row["level"] == "network"
             and row["condition"] == "open" and row["estimator"] == "aec"
             and row["band"] == "beta" and row["network"] == "Vis"
             and row["quantity"] == "auc_degree" and row["aqvar"] == "detail"]
    require(len(focal) == 1 and round(focal[0]["r"], 3) == .539, "Focal visual result differs")
    with (DATA / "results_regression_summary.csv").open(encoding="utf-8-sig", newline="") as stream:
        regression = list(csv.DictReader(stream))
    units = [row for row in regression if row["family"] == "confirm"]
    require(len(regression) == 8 and len(units) == 7, "Unexpected follow-up unit counts")
    for row in units:
        require(abs(float(row["q2_unit"]) - float(row["q2_base"]) - float(row["dq2"])) < 1e-12,
                "Conditional prediction difference is inconsistent")
    return {
        "status": "PASS",
        "scope": "Aggregate-result audit; participant-level estimation is not rerun",
        "tests": len(rows), "local_families": families,
        "local_hits": len(hits), "hit_families": hit_families,
        "hits_by_AQ_score": dict(Counter(row["aqvar"] for row in hits)),
        "maximum_local_q_difference": error,
        "band_pooled_families": pooled_families, "band_pooled_survivors": 0,
        "minimum_band_pooled_q": min(pooled), "studywide_survivors": 0,
        "minimum_studywide_q": min(studywide),
        "focal_visual_result": {key: focal[0][key] for key in ("r", "ci_lo", "ci_hi", "p", "q4", "n")},
        "reported_followup_units": len(units),
        "historical_extra_regression_rows": len(regression) - len(units),
    }


if __name__ == "__main__":
    try:
        print(json.dumps(audit(), indent=2))
    except (ValueError, KeyError, OSError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
