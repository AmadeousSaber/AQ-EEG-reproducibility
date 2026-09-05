"""Standard-library tests for the public aggregate audit."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from verify_public_results import audit, bh_adjust, family_key


class PublicVerificationTests(unittest.TestCase):
    def test_bh_known_values_and_ties(self):
        expected = [.04, .04, .04, .9]
        for actual, target in zip(bh_adjust([.01, .02, .03, .9]), expected, strict=True):
            self.assertAlmostEqual(actual, target)
        self.assertEqual(bh_adjust([0, 0, 1]), [0, 0, 1])
        self.assertEqual(bh_adjust([]), [])

    def test_bh_rejects_invalid_values(self):
        for value in (float("nan"), float("inf"), -.1, 1.1):
            with self.assertRaises(ValueError):
                bh_adjust([value])

    def test_band_pooling_retains_condition_estimator_and_scale(self):
        row = dict(condition="open", estimator="aec", level="global", band="beta",
                   metric="aec_beta_auc_gcc")
        self.assertEqual(family_key(row), ("open", "aec", "global", "auc_gcc", "beta"))
        self.assertEqual(family_key(row, True), ("open", "aec", "global", "auc_gcc"))

    def test_frozen_public_results(self):
        result = audit()
        self.assertEqual(result["tests"], 20928)
        self.assertEqual(result["local_hits"], 14)
        self.assertEqual(result["reported_followup_units"], 7)


if __name__ == "__main__":
    unittest.main()
