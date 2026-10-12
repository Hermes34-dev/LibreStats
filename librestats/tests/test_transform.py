"""Check librestats.transform against scipy.

Run: ../.venv-dev/bin/python -m unittest tests.test_transform   (from librestats/)
"""

import math
import os
import sys
import unittest

import numpy as np
from scipy import special, stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import transform as t  # noqa: E402

RNG = np.random.default_rng(2026)
SKEWED = RNG.lognormal(1.0, 0.8, 60).tolist()
WITH_NEG = (RNG.standard_normal(60) * 3 + RNG.exponential(4, 60)).tolist()


def p(res, label):
    return dict(res.params)[label]


class TestLogSqrtLogit(unittest.TestCase):
    def test_log_offsets(self):
        vals = [0.0, 2.0, None, 8.0]
        self.assertEqual(t.log_transform(vals, "e", "one").values,
                         [0.0, math.log(3), None, math.log(9)])
        self.assertAlmostEqual(t.log_transform(vals, "10", "half_min").values[1], math.log10(3))
        with self.assertRaises(t.TransformError):
            t.log_transform(vals)  # log(0)

    def test_sqrt(self):
        self.assertEqual(t.sqrt_transform([4.0, None, 0.0]).values, [2.0, None, 0.0])
        with self.assertRaises(t.TransformError):
            t.sqrt_transform([-1.0, 4.0])

    def test_logit(self):
        x = [0.1, 0.5, 0.9, None]
        np.testing.assert_allclose(t.logit_transform(x).values[:3], special.logit(x[:3]))
        pct = t.logit_transform([10.0, 50.0, 90.0], percent=True).values
        np.testing.assert_allclose(pct, special.logit([0.1, 0.5, 0.9]))

    def test_logit_squeeze(self):
        x = [0.0, 0.3, 1.0, 0.6]
        n = 4
        expected = special.logit([(v * (n - 1) + 0.5) / n for v in x])
        res = t.logit_transform(x)
        np.testing.assert_allclose(res.values, expected)
        self.assertTrue(res.notes)
        with self.assertRaises(t.TransformError):
            t.logit_transform(x, squeeze=False)


class TestPowerTransforms(unittest.TestCase):
    def test_boxcox_matches_scipy(self):
        y_ref, lam_ref, (lo_ref, hi_ref) = stats.boxcox(SKEWED, alpha=0.05)
        res = t.boxcox(SKEWED + [None])
        self.assertAlmostEqual(p(res, "Lambda (ML)"), lam_ref, places=6)
        self.assertAlmostEqual(p(res, "95% CI lower"), lo_ref, places=6)
        self.assertAlmostEqual(p(res, "95% CI upper"), hi_ref, places=6)
        np.testing.assert_allclose(res.values[:-1], y_ref, rtol=1e-6)
        self.assertIsNone(res.values[-1])

    def test_boxcox_llf_matches_scipy(self):
        for lam in (-1.5, 0.0, 0.3, 2.0):
            self.assertAlmostEqual(t.boxcox_llf(lam, SKEWED), stats.boxcox_llf(lam, SKEWED),
                                   places=8)

    def test_boxcox_fixed_lambda_and_errors(self):
        res = t.boxcox(SKEWED, lmbda=0)
        np.testing.assert_allclose(res.values, np.log(SKEWED))
        with self.assertRaises(t.TransformError):
            t.boxcox([-1.0, 2.0, 3.0])
        shifted = t.boxcox([0.0, 1.0, 2.0, 5.0, 9.0], offset="one")
        self.assertEqual(p(shifted, "Offset added"), 1.0)

    def test_yeojohnson_matches_scipy(self):
        y_ref, lam_ref = stats.yeojohnson(WITH_NEG)
        res = t.yeojohnson(WITH_NEG)
        self.assertAlmostEqual(p(res, "Lambda (ML)"), lam_ref, places=5)
        np.testing.assert_allclose(res.values, stats.yeojohnson(WITH_NEG, lmbda=p(res, "Lambda (ML)")),
                                   rtol=1e-10)
        for lam in (-0.5, 0.0, 1.0, 2.0, 2.7):
            self.assertAlmostEqual(t.yeojohnson_llf(lam, WITH_NEG),
                                   stats.yeojohnson_llf(lam, WITH_NEG), places=8)
            np.testing.assert_allclose(t.yeojohnson(WITH_NEG, lmbda=lam).values,
                                       stats.yeojohnson(WITH_NEG, lmbda=lam), rtol=1e-10, atol=1e-12)


class TestScaling(unittest.TestCase):
    def test_methods(self):
        x = np.array(SKEWED)
        mean, sd = x.mean(), x.std(ddof=1)
        vals = SKEWED + [None]
        np.testing.assert_allclose(t.scale(vals, "centre").values[:-1], x - mean)
        np.testing.assert_allclose(t.scale(vals, "z").values[:-1], stats.zscore(x, ddof=1))
        np.testing.assert_allclose(t.scale(vals, "range").values[:-1],
                                   (x - x.min()) / (x.max() - x.min()))
        np.testing.assert_allclose(t.scale(vals, "pareto").values[:-1], (x - mean) / np.sqrt(sd))
        self.assertIsNone(t.scale(vals, "z").values[-1])

    def test_constant(self):
        with self.assertRaises(t.TransformError):
            t.scale([2.0, 2.0, 2.0], "z")


class TestRanksAndNormalScores(unittest.TestCase):
    def test_ranks_match_scipy(self):
        x = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0]
        for ties in t.TIES:
            self.assertEqual(t.ranks(x, ties).values, stats.rankdata(x, method=ties).tolist())
        self.assertEqual(t.ranks([2.0, None, 1.0]).values, [2.0, None, 1.0])

    def test_normal_scores(self):
        x = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, None]
        r = stats.rankdata(x[:-1])
        n = len(r)
        for method, a in t.NORMAL_SCORES.items():
            ref = stats.norm.ppf((r - a) / (n + 1 - 2 * a))
            res = t.normal_scores(x, method).values
            np.testing.assert_allclose(res[:-1], ref, rtol=1e-12)
            self.assertIsNone(res[-1])


if __name__ == "__main__":
    unittest.main()
