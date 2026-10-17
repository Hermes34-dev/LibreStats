"""Check librestats.assumptions (Shapiro-Wilk, Levene) against scipy and by hand."""

import math
import os
import sys
import unittest

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import assumptions as a  # noqa: E402

RNG = np.random.default_rng(31)
SAMPLES = {
    "normal": lambda n: RNG.normal(10, 2, n),
    "exponential": lambda n: RNG.exponential(1, n),
    "uniform": lambda n: RNG.uniform(0, 1, n),
    "t3": lambda n: RNG.standard_t(3, n),
}


class TestShapiroWilk(unittest.TestCase):
    def test_matches_scipy(self):
        for n in (3, 4, 5, 6, 7, 11, 12, 20, 50, 200, 1000, 5000):
            for dist, draw in SAMPLES.items():
                x = draw(n)
                res, ref = a.shapiro_wilk(list(x)), stats.shapiro(x)
                self.assertAlmostEqual(res.statistic, ref.statistic, delta=1e-10, msg=(n, dist))
                self.assertAlmostEqual(res.p, ref.pvalue, delta=1e-9, msg=(n, dist))
                self.assertEqual(res.n, n)

    def test_coefficients(self):
        for n in (3, 4, 5, 6, 10, 51, 1000):
            c = a.shapiro_coefficients(n)
            self.assertEqual(len(c), n)
            self.assertAlmostEqual(sum(v * v for v in c), 1.0, delta=1e-12)
            self.assertTrue(all(abs(c[i] + c[n - 1 - i]) < 1e-15 for i in range(n)))

    def test_n3_exact(self):
        # equally spaced values: W = 1, p = 1; one value far off: W near 0.75
        res = a.shapiro_wilk([1.0, 2.0, 3.0])
        self.assertAlmostEqual(res.statistic, 1.0, delta=1e-12)
        self.assertAlmostEqual(res.p, 1.0, delta=1e-12)
        res = a.shapiro_wilk([1.0, 1.001, 100.0])
        ref = stats.shapiro([1.0, 1.001, 100.0])
        self.assertAlmostEqual(res.p, ref.pvalue, delta=1e-9)

    def test_rejects_skewed_large_sample(self):
        self.assertLess(a.shapiro_wilk(list(RNG.exponential(1, 300))).p, 1e-6)

    def test_order_and_missing_ignored(self):
        x = list(RNG.normal(0, 1, 25))
        res = a.shapiro_wilk(x[::-1] + [None])
        self.assertAlmostEqual(res.statistic, stats.shapiro(x).statistic, delta=1e-12)

    def test_errors_and_notes(self):
        for x in ([], [1.0], [1.0, 2.0], [4.0] * 10):
            with self.assertRaises(a.AssumptionError, msg=x):
                a.shapiro_wilk(x)
        self.assertTrue(any("5000" in n for n in a.shapiro_wilk(list(RNG.normal(0, 1, 5001)))
                            .notes))
        self.assertTrue(any("Tied" in n for n in a.shapiro_wilk([1.0, 2.0, 2.0, 3.0, 5.0])
                            .notes))


class TestLevene(unittest.TestCase):
    def test_by_hand(self):
        # median-centred |deviations|: [1, 0, 1] and [2, 0, 2]; F = 4 * (2/3) / (10/3) = 0.8
        res = a.levene([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0]])
        self.assertAlmostEqual(res.statistic, 0.8, delta=1e-12)
        self.assertEqual(dict(res.params), {"df1": 1, "df2": 4})
        self.assertAlmostEqual(res.p, stats.f.sf(0.8, 1, 4), delta=1e-12)

    def test_matches_scipy(self):
        for center in ("median", "mean"):
            for sizes, sds in (((10, 12), (1, 1)), ((8, 15, 9), (1, 2, 1.5)),
                               ((30, 31, 40, 25, 7), (1, 1, 1.2, 0.8, 3)), ((3, 3), (1, 5))):
                groups = [RNG.normal(5, s, m) for m, s in zip(sizes, sds)]
                res = a.levene([list(g) for g in groups], center)
                ref = stats.levene(*groups, center=center)
                self.assertAlmostEqual(res.statistic, ref.statistic, delta=1e-10,
                                       msg=(center, sizes))
                self.assertAlmostEqual(res.p, ref.pvalue, delta=1e-10, msg=(center, sizes))

    def test_skewed_groups(self):
        groups = [RNG.exponential(s, 40) for s in (1, 1, 3)]
        for center in ("median", "mean"):
            ref = stats.levene(*groups, center=center)
            self.assertAlmostEqual(a.levene([list(g) for g in groups], center).p, ref.pvalue,
                                   delta=1e-10)

    def test_errors(self):
        for groups in ([[1.0, 2.0, 3.0]], [[1.0, 2.0], [3.0]], [[1.0, 3.0], [2.0, 6.0]]):
            with self.assertRaises(a.AssumptionError, msg=groups):
                a.levene(groups)
        with self.assertRaises(a.AssumptionError):
            a.levene([[1.0, 2.0, 3.0], [1.0, 5.0, 9.0]], center="trimmed")


if __name__ == "__main__":
    unittest.main()
