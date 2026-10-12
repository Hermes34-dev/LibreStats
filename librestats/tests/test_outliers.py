"""Check librestats.outliers against numpy/scipy, published values and simulation.

Run: /path/to/.venv-dev/bin/python -m unittest discover -s tests   (from librestats/)
"""

import os
import sys
import unittest

import numpy as np
from scipy import stats
from scipy.spatial.distance import mahalanobis as scipy_mahalanobis

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import outliers as o  # noqa: E402

RNG = np.random.default_rng(12345)
DATA = [2.1, 2.4, 2.2, 2.3, 2.5, 2.2, 2.4, 9.8, 2.3, 2.6, 2.1, -3.0, 2.4]


class TestTukey(unittest.TestCase):
    def test_fences_match_numpy_percentile(self):
        res = o.tukey_fences(DATA, k=1.5)
        q1, q3 = np.percentile(DATA, [25, 75])  # numpy default = type 7
        summary = dict(res.summary)
        self.assertAlmostEqual(summary["Q1"], q1, places=12)
        self.assertAlmostEqual(summary["Q3"], q3, places=12)
        self.assertAlmostEqual(summary["Lower fence"], q1 - 1.5 * (q3 - q1), places=12)
        self.assertEqual(sorted(f.index for f in res.flags), [7, 11])
        self.assertEqual({f.index: f.side for f in res.flags}, {7: "high", 11: "low"})

    def test_k3_is_stricter(self):
        self.assertLessEqual(len(o.tukey_fences(DATA, 3.0).flags),
                             len(o.tukey_fences(DATA, 1.5).flags))


class TestMAD(unittest.TestCase):
    def test_matches_scipy(self):
        res = o.mad_zscores(DATA, 2.5)
        mad = stats.median_abs_deviation(DATA, scale="normal")  # 1/Phi^-1(0.75) = 1.482602...
        self.assertAlmostEqual(dict(res.summary)["MAD (scaled, b = 1.4826)"] / mad, 1, places=5)
        med = np.median(DATA)
        scaled = 1.4826 * stats.median_abs_deviation(DATA)
        expected = [i for i, x in enumerate(DATA) if abs(x - med) / scaled > 2.5]
        self.assertEqual(sorted(f.index for f in res.flags), expected)

    def test_zero_mad_raises(self):
        with self.assertRaises(o.OutlierError):
            o.mad_zscores([1, 1, 1, 1, 5])


class TestGrubbs(unittest.TestCase):
    def test_published_critical_values(self):
        # Two-sided alpha = 0.05 (e.g. ASTM E178 / NIST handbook tables)
        for n, g in ((3, 1.1543), (10, 2.2900), (20, 2.7082)):
            self.assertAlmostEqual(o.grubbs_critical(n, 0.05), g, places=3)

    def test_critical_matches_scipy_formula(self):
        for n in (3, 7, 25, 100):
            for a in (0.1, 0.05, 0.01):
                t = stats.t.isf(a / (2 * n), n - 2)
                ref = (n - 1) / np.sqrt(n) * np.sqrt(t * t / (n - 2 + t * t))
                self.assertAlmostEqual(o.grubbs_critical(n, a), ref, places=9)

    def test_pvalue_consistent_with_critical(self):
        for n in (5, 15, 40):
            g = o.grubbs_critical(n, 0.05)
            self.assertAlmostEqual(o.grubbs_pvalue(g, n), 0.05, places=9)

    def test_detects_outlier_and_iterates(self):
        single = o.grubbs(DATA, 0.05)
        self.assertEqual([f.index for f in single.flags], [7])
        repeated = o.grubbs(DATA, 0.05, iterate=True)
        self.assertEqual([f.index for f in repeated.flags], [7, 11])
        self.assertFalse(repeated.steps[-1]["outlier"])

    def test_null_rejection_rate(self):
        n, reps = 12, 40000
        x = RNG.standard_normal((reps, n))
        g = np.max(np.abs(x - x.mean(1, keepdims=True)), 1) / x.std(1, ddof=1)
        rate = np.mean(g > o.grubbs_critical(n, 0.05))
        self.assertAlmostEqual(rate, 0.05, delta=0.005)


class TestDixon(unittest.TestCase):
    def test_classic_q_test_example(self):
        # Ten replicate measurements; 0.167 is suspicious but not significant.
        data = [0.189, 0.167, 0.187, 0.183, 0.186, 0.182, 0.181, 0.184, 0.181, 0.177]
        res = o.dixon(data, 0.05)
        step = res.steps[0]
        self.assertEqual(step["ratio"], "r11")
        self.assertAlmostEqual(step["r"], (0.177 - 0.167) / (0.187 - 0.167), places=12)
        self.assertEqual(step["index"], 1)
        self.assertFalse(res.flags)

    def test_obvious_outlier(self):
        res = o.dixon([10.1, 10.3, 10.2, 10.4, 10.2, 14.9], 0.05)
        self.assertEqual([(f.index, f.side) for f in res.flags], [(5, "high")])
        self.assertEqual(res.steps[0]["p"], "p < 0.01")

    def test_null_rejection_rate(self):
        # Pure-Python test on normal data should reject ~alpha of the time.
        for n in (5, 9, 12, 25):
            reps = 6000
            rejections = sum(bool(o.dixon(list(row), 0.05).flags)
                             for row in RNG.standard_normal((reps, n)))
            self.assertAlmostEqual(rejections / reps, 0.05, delta=0.012, msg=f"n={n}")

    def test_size_limits(self):
        with self.assertRaises(o.OutlierError):
            o.dixon(list(range(31)))
        with self.assertRaises(o.OutlierError):
            o.dixon([1.0, 2.0])


class TestMahalanobis(unittest.TestCase):
    def test_matches_scipy(self):
        x = RNG.multivariate_normal([0, 5, -2], [[1, .6, .2], [.6, 2, .3], [.2, .3, .5]], 60)
        x[10] = [3.0, 1.0, 0.5]  # off the correlation structure
        rows = x.tolist()
        res = o.mahalanobis(rows, alpha=0.001)
        vi = np.linalg.inv(np.cov(x, rowvar=False))
        mu = x.mean(0)
        for i, (d2, p) in enumerate(res.scores):
            ref = scipy_mahalanobis(x[i], mu, vi) ** 2
            self.assertAlmostEqual(d2, ref, places=8)
            self.assertAlmostEqual(p, stats.chi2.sf(ref, 3), places=10)
        self.assertAlmostEqual(dict(res.summary)["Critical D²"], stats.chi2.isf(0.001, 3), places=8)
        self.assertIn(10, [f.index for f in res.flags])

    def test_singular_raises(self):
        rows = [[i, 2 * i + 1] for i in range(10)]
        with self.assertRaises(o.OutlierError):
            o.mahalanobis(rows)

    def test_unreachable_cutoff_is_noted(self):
        x = RNG.standard_normal((13, 3)).tolist()
        self.assertTrue(o.mahalanobis(x, 0.001).notes)
        big = RNG.standard_normal((200, 3)).tolist()
        self.assertFalse(o.mahalanobis(big, 0.001).notes)

    def test_too_few_rows(self):
        with self.assertRaises(o.OutlierError):
            o.mahalanobis([[1, 2, 3], [2, 3, 1], [3, 1, 2], [1, 1, 1]])


if __name__ == "__main__":
    unittest.main()
