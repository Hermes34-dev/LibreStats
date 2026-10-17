"""Check librestats.distributions against scipy (dev-only reference).

Run: ../.venv-dev/bin/python -m unittest discover -s tests   (from librestats/)
"""

import os
import sys
import unittest

from scipy import special, stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import distributions as d  # noqa: E402

DFS = (1, 2, 3, 5, 8, 13, 30, 100, 1000)
PS = (0.5, 0.2, 0.05, 0.01, 1e-3, 1e-5, 1e-8)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


class TestSpecialFunctions(unittest.TestCase):
    def test_betainc(self):
        for a in (0.5, 1, 2.5, 10, 50):
            for b in (0.5, 1, 3, 20):
                for x in (1e-6, 0.01, 0.3, 0.5, 0.77, 0.999):
                    self.assertLess(rel(d.betainc(a, b, x), special.betainc(a, b, x)), 1e-10,
                                    (a, b, x))

    def test_gammainc_upper(self):
        for a in (0.5, 1, 2.5, 10, 100):
            for x in (1e-4, 0.5, 3, 10, 50, 150):
                ref = special.gammaincc(a, x)
                if ref < 1e-290:
                    continue
                self.assertLess(rel(d.gammainc_upper(a, x), ref), 1e-10, (a, x))


class TestT(unittest.TestCase):
    def test_sf(self):
        for df in DFS:
            for t in (-3.0, -0.5, 0.0, 0.7, 2.0, 5.0, 20.0):
                self.assertLess(rel(d.t_sf(t, df), stats.t.sf(t, df)), 1e-9, (t, df))

    def test_isf(self):
        for df in DFS:
            for p in PS + (0.9, 0.99):
                self.assertLess(rel(d.t_isf(p, df), stats.t.isf(p, df)) if p != 0.5
                                else abs(d.t_isf(p, df)), 1e-9, (p, df))


class TestChi2(unittest.TestCase):
    def test_sf(self):
        for df in DFS:
            for x in (0.01, 0.5, 1, 5, 20, 80):
                ref = stats.chi2.sf(x, df)
                if ref < 1e-290:
                    continue
                self.assertLess(rel(d.chi2_sf(x, df), ref), 1e-9, (x, df))

    def test_isf(self):
        for df in DFS:
            for p in PS + (0.9, 0.99):
                self.assertLess(rel(d.chi2_isf(p, df), stats.chi2.isf(p, df)), 1e-9, (p, df))


class TestF(unittest.TestCase):
    def test_sf(self):
        for df1 in (1, 2, 4, 9, 30):
            for df2 in (1, 3, 10, 57, 500):
                for x in (0.0, 0.05, 0.5, 1.0, 2.7, 8.0, 40.0):
                    ref = stats.f.sf(x, df1, df2)
                    self.assertLess(rel(d.f_sf(x, df1, df2), ref), 1e-9, (x, df1, df2))


class TestNormal(unittest.TestCase):
    def test_cdf_and_ppf(self):
        for x in (-30, -8, -2.5, -0.3, 0, 0.7, 3, 8):
            self.assertLess(rel(d.norm_cdf(x), stats.norm.cdf(x)), 1e-12, x)
        for p in (1e-300, 1e-15, 1e-6, 0.01, 0.2, 0.5, 0.8, 0.975, 1 - 1e-10):
            ref = stats.norm.ppf(p)
            self.assertLess(abs(d.norm_ppf(p) - ref), 1e-12 * max(1, abs(ref)), p)


if __name__ == "__main__":
    unittest.main()
