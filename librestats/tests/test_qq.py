"""Check librestats.qq against scipy/numpy and R's qqnorm / qqline formulas."""

import os
import sys
import unittest

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import qq  # noqa: E402

RNG = np.random.default_rng(5)


class TestPositions(unittest.TestCase):
    def test_r_ppoints_n5(self):
        # R: qnorm(ppoints(5)) = -1.1797611 -0.4972006 0.0000000 0.4972006 1.1797611
        z = qq.plotting_positions(5, "blom")
        for got, ref in zip(z, (-1.1797611, -0.4972006, 0.0, 0.4972006, 1.1797611)):
            self.assertAlmostEqual(got, ref, delta=1e-7)

    def test_all_methods_match_scipy(self):
        for method, a in (("blom", 3 / 8), ("vdw", 0.0), ("tukey", 1 / 3), ("rankit", 0.5)):
            for n in (3, 10, 11, 250):
                ref = stats.norm.ppf((np.arange(1, n + 1) - a) / (n + 1 - 2 * a))
                np.testing.assert_allclose(qq.plotting_positions(n, method), ref, atol=1e-12)


class TestQQ(unittest.TestCase):
    def test_sorted_sample_and_ties(self):
        res = qq.qq_normal([3.0, 1.0, 2.0, 2.0, None], line="none")
        self.assertEqual(res.sample, [1.0, 2.0, 2.0, 3.0])
        self.assertEqual(len(set(res.theoretical)), 4)   # ties get their own positions
        self.assertIsNone(res.slope)
        self.assertEqual(res.line(), [])

    def test_quartile_line_as_r_qqline(self):
        x = RNG.gamma(2, 3, 57)
        res = qq.qq_normal(list(x))
        q1, q3 = np.percentile(x, [25, 75])          # type 7, as R's quantile default
        z1, z3 = stats.norm.ppf([0.25, 0.75])
        slope = (q3 - q1) / (z3 - z1)
        self.assertAlmostEqual(res.slope, slope, delta=1e-12)
        self.assertAlmostEqual(res.intercept, q1 - slope * z1, delta=1e-12)
        self.assertAlmostEqual(res.line()[0], res.intercept + res.slope * res.theoretical[0],
                               delta=1e-12)

    def test_mean_sd_line(self):
        x = RNG.normal(4, 2, 40)
        res = qq.qq_normal(list(x), line="meansd")
        self.assertAlmostEqual(res.intercept, x.mean(), delta=1e-12)
        self.assertAlmostEqual(res.slope, x.std(ddof=1), delta=1e-12)

    def test_normal_data_lie_near_line(self):
        x = RNG.normal(10, 2, 2000)
        res = qq.qq_normal(list(x))
        self.assertAlmostEqual(res.intercept, 10, delta=0.2)
        self.assertAlmostEqual(res.slope, 2, delta=0.2)

    def test_errors_and_notes(self):
        for x in ([], [1.0, 2.0]):
            with self.assertRaises(qq.QQError):
                qq.qq_normal(x)
        with self.assertRaises(qq.QQError):
            qq.qq_normal([1.0, 2.0, 3.0], method="filliben")
        with self.assertRaises(qq.QQError):
            qq.qq_normal([1.0, 2.0, 3.0], line="loess")
        res = qq.qq_normal([5.0] * 4)
        self.assertEqual(res.slope, 0)
        self.assertTrue(res.notes)


if __name__ == "__main__":
    unittest.main()
