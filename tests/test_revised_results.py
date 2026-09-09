"""Behavioral tests for version 1.1.0 aggregate verification."""
from copy import deepcopy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from verify_revised_results import EXPECTED_ROWS, PREFIX, ROOT, STABILITY, audit, check_tables, plus_one, read_tables, wilson


class RevisedVerificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables = read_tables()

    def test_complete_snapshot(self):
        report = audit()
        self.assertEqual(report["revised_csv_files"], 8)
        self.assertEqual(report["permutation_effect_statistic_rows"], 28)

    def test_plus_one_boundaries(self):
        self.assertEqual(plus_one(0, 999), .001)
        self.assertEqual(plus_one(999, 999), 1)
        for count, trials in ((-1, 999), (1000, 999), (0, 0), (1.5, 999)):
            with self.assertRaises(ValueError):
                plus_one(count, trials)

    def test_wilson_zero_and_full(self):
        lo, hi = wilson(0, 1000)
        self.assertAlmostEqual(lo, 0)
        self.assertGreater(hi, 0)
        lo, hi = wilson(1000, 1000)
        self.assertLess(lo, 1)
        self.assertAlmostEqual(hi, 1)

    def test_wrong_permutation_count_is_rejected(self):
        tables = deepcopy(self.tables)
        tables[PREFIX + "results_regression_v2_permutations.csv"][0]["exceedances"] = "9000"
        with self.assertRaises(ValueError):
            check_tables(tables)

    def test_corrupt_bh_is_rejected(self):
        tables = deepcopy(self.tables)
        tables[PREFIX + "results_regression_v2_permutations.csv"][0]["q_BH_seven"] = ".8"
        with self.assertRaises(ValueError):
            check_tables(tables)

    def test_interval_and_summary_mismatch_are_rejected(self):
        for column, value in (("dR2_boot_lo", "0.8"), ("dq2", "4"), ("family", "confirm")):
            tables = deepcopy(self.tables)
            tables[PREFIX + "results_regression_v2_summary.csv"][0][column] = value
            with self.assertRaises(ValueError):
                check_tables(tables)

    def test_inconsistent_simulation_rate_is_rejected(self):
        tables = deepcopy(self.tables)
        tables[STABILITY + "interaction_null_summary.csv"][0]["rejection_rate"] = ".5"
        with self.assertRaises(ValueError):
            check_tables(tables)

    def test_duplicate_row_is_rejected(self):
        tables = deepcopy(self.tables)
        rows = tables[PREFIX + "results_regression_v2_permutations.csv"]
        rows[1] = dict(rows[0])
        with self.assertRaises(ValueError):
            check_tables(tables)

    def test_changed_file_hash_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            shutil.copyfile(ROOT / "revised_manifest.json", root / "revised_manifest.json")
            for relative in EXPECTED_ROWS:
                destination = root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT / relative, destination)
            changed = root / (PREFIX + "results_regression_v2_summary.csv")
            content = changed.read_bytes()
            changed.write_bytes(b"X" + content[1:])
            with self.assertRaisesRegex(ValueError, "Checksum mismatch"):
                read_tables(root)

    def test_manifest_path_traversal_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = json.loads((ROOT / "revised_manifest.json").read_text(encoding="utf-8"))
            manifest["files"][0]["path"] = "../private.csv"
            (root / "revised_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unexpected revised manifest paths"):
                read_tables(root)

    def test_nonfinite_p_is_rejected(self):
        tables = deepcopy(self.tables)
        tables[PREFIX + "results_regression_v2_permutations.csv"][0]["p_plus_one"] = "NaN"
        with self.assertRaisesRegex(ValueError, "Nonfinite"):
            check_tables(tables)


if __name__ == "__main__":
    unittest.main()
