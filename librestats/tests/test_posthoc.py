"""Check the studentized range and Dunnett distributions and librestats.posthoc
(Tukey HSD, Games-Howell, Dunnett, Dunn) against scipy, statsmodels and numpy."""

import math
import os
import sys
import unittest

import numpy as np
from scipy import special, stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import distributions as D  # noqa: E402
from librestats import posthoc as ph  # noqa: E402

RNG = np.random.default_rng(59)


def rel(a, b):
    return abs(a - b) / abs(b)


def brute_tukey_sf(q, k, df, ns=1501, nz=6001):
    """Independent check: trapezoid rule in s and z on fine grids (numpy)."""
    s = np.linspace(1e-6, 2.5, ns)
    logf = (np.log(2) + df / 2 * np.log(df / 2) - special.gammaln(df / 2) + (df - 1) * np.log(s)
            - df * s * s / 2)
    inner = []
    for w in q * s:
        z = np.linspace(-w / 2 - 12, 12, nz)
        a, c = special.ndtr(-z), special.ndtr(-(z + w))
        bracket = np.where(c >= a, 1.0, -np.expm1((k - 1) * np.log1p(-c / a)))
        inner.append(np.trapezoid(k * stats.norm.pdf(z) * a ** (k - 1) * bracket, z))
    return np.trapezoid(np.exp(logf) * np.array(inner), s)


class TestStudentizedRange(unittest.TestCase):
    def test_two_means_is_t(self):
        # Q / sqrt(2) = |T| for k = 2: exact down to tiny p-values
        for df in (1, 1.7, 3, 10, 25.5, 200, 1000, 1e5):
            for q in (0.3, 1, 2, 4, 8, 12, 20):
                ref = 2 * stats.t.sf(q / math.sqrt(2), df)
                self.assertLess(rel(D.tukey_sf(q, 2, df), ref), 1e-9, msg=(df, q))

    def test_matches_scipy(self):
        # scipy's studentized_range is accurate to ~1e-8 only for moderate p and df
        for k in (3, 5, 10, 20):
            for df in (1, 1.7, 3, 10, 25.5, 60, 200):
                for q in (0.3, 1, 2, 3, 4, 6, 8):
                    ref = stats.studentized_range.sf(q, k, df)
                    if ref > 1e-6:
                        self.assertLess(rel(D.tukey_sf(q, k, df), ref), 1e-7, msg=(k, df, q))

    def test_far_tail_brute_force(self):
        for k, df, q in ((3, 60, 12), (10, 30, 9)):
            self.assertLess(rel(D.tukey_sf(q, k, df), brute_tukey_sf(q, k, df)), 1e-7)

    def test_quantiles(self):
        for k, df, p in ((3, 10, .05), (5, 20, .05), (10, 60, .01), (4, 2.5, .05), (2, 7, .1)):
            q = D.tukey_isf(p, k, df)
            self.assertLess(rel(q, stats.studentized_range.ppf(1 - p, k, df)), 1e-9)
        # published table: q(0.05; 3, 10) = 3.877, q(0.05; 5, 20) = 4.232
        self.assertAlmostEqual(D.tukey_isf(0.05, 3, 10), 3.877, delta=5e-4)
        self.assertAlmostEqual(D.tukey_isf(0.05, 5, 20), 4.232, delta=5e-4)


class TestDunnettDistribution(unittest.TestCase):
    def test_one_treatment_is_t(self):
        lam = [math.sqrt(8 / 18)]
        for df in (3, 17, 100):
            for t in (0.5, 2, 5):
                self.assertLess(rel(D.dunnett_sf(t, df, lam), 2 * stats.t.sf(t, df)), 1e-10)
                self.assertLess(rel(D.dunnett_sf(t, df, lam, False), stats.t.sf(t, df)), 1e-10)

    def test_matches_multivariate_t(self):
        rng = np.random.default_rng(1)
        for ns, n0, df in (((10, 10, 10), 10, 27), ((5, 12, 8, 20), 9, 50)):
            lam = [math.sqrt(n / (n + n0)) for n in ns]
            corr = np.outer(lam, lam)
            np.fill_diagonal(corr, 1)
            mvt = stats.multivariate_t(np.zeros(len(ns)), corr, df=df)
            for t in (1.5, 2.4, 3.2):
                upper = np.full(len(ns), t)
                two = 1 - mvt.cdf(upper, lower_limit=-upper, maxpts=10 ** 6, random_state=rng)
                one = 1 - mvt.cdf(upper, maxpts=10 ** 6, random_state=rng)
                self.assertLess(rel(D.dunnett_sf(t, df, lam), two), 1e-4, msg=(ns, t))
                self.assertLess(rel(D.dunnett_sf(t, df, lam, False), one), 1e-4, msg=(ns, t))

    def test_published_critical_values(self):
        # Dunnett (1964) tables, two-sided 0.05, equal n: 2 treatments df 20: 2.38;
        # 4 treatments df 30: 2.58; one-sided 0.05, 3 treatments df 24: 2.17
        half = [math.sqrt(0.5)]
        self.assertAlmostEqual(D.dunnett_isf(0.05, 20, half * 2), 2.38, delta=0.006)
        self.assertAlmostEqual(D.dunnett_isf(0.05, 30, half * 4), 2.58, delta=0.006)
        self.assertAlmostEqual(D.dunnett_isf(0.05, 24, half * 3, False), 2.17, delta=0.006)
        c = D.dunnett_isf(0.05, 27, half * 3)
        self.assertAlmostEqual(D.dunnett_sf(c, 27, half * 3), 0.05, delta=1e-10)


def make_groups(sizes, sds, shift=0.4):
    return [RNG.normal(10 + i * shift, s, n) for i, (n, s) in enumerate(zip(sizes, sds))]


class TestTukeyHSD(unittest.TestCase):
    def test_matches_scipy(self):
        for sizes in ((6, 6, 6), (5, 9, 14, 7), (10, 12)):
            groups = make_groups(sizes, [1] * len(sizes))
            res = ph.tukey_hsd([list(g) for g in groups], 0.05)
            ref = stats.tukey_hsd(*groups)
            ci = ref.confidence_interval(0.95)
            for c in res.comparisons:
                i, j = c.first, c.second
                self.assertAlmostEqual(c.estimate, ref.statistic[i, j], delta=1e-12)
                self.assertLess(rel(c.p, ref.pvalue[i, j]), 1e-6, msg=(sizes, i, j))
                self.assertAlmostEqual(c.ci[0], ci.low[i, j], delta=1e-8)
                self.assertAlmostEqual(c.ci[1], ci.high[i, j], delta=1e-8)
            self.assertEqual(len(res.comparisons), len(sizes) * (len(sizes) - 1) // 2)


class TestGamesHowell(unittest.TestCase):
    def test_formula(self):
        groups = make_groups((5, 12, 8, 20), (0.5, 2, 1, 3))
        res = ph.games_howell([list(g) for g in groups], 0.05)
        k = len(groups)
        for c in res.comparisons:
            a, b = groups[c.first], groups[c.second]
            va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
            df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
            t = (a.mean() - b.mean()) / math.sqrt(va + vb)
            self.assertAlmostEqual(c.statistic, t, delta=1e-10)
            self.assertAlmostEqual(c.df, df, delta=1e-9)
            p = stats.studentized_range.sf(abs(t) * math.sqrt(2), k, df)
            self.assertLess(rel(c.p, p), 1e-6)
            half = stats.studentized_range.ppf(0.95, k, df) / math.sqrt(2) * math.sqrt(va + vb)
            self.assertAlmostEqual(c.ci[1] - c.ci[0], 2 * half, delta=1e-7)


class TestDunnett(unittest.TestCase):
    def test_matches_scipy(self):
        # scipy.stats.dunnett integrates by randomized quasi-Monte Carlo (~1e-4)
        for sizes, control in (((8, 8, 8, 8), 0), ((10, 5, 12, 7), 2)):
            groups = make_groups(sizes, [1] * len(sizes), 0.6)
            for alt in ("two-sided", "greater", "less"):
                res = ph.dunnett([list(g) for g in groups], control, alt)
                others = [g for i, g in enumerate(groups) if i != control]
                ref = stats.dunnett(*others, control=groups[control], alternative=alt,
                                    random_state=np.random.default_rng(3))
                ci = ref.confidence_interval(0.95)
                for m, c in enumerate(res.comparisons):
                    self.assertAlmostEqual(c.statistic, ref.statistic[m], delta=1e-10)
                    self.assertAlmostEqual(c.p, ref.pvalue[m], delta=2e-4, msg=(sizes, alt, m))
                    for got, want in zip(c.ci, (ci.low[m], ci.high[m])):
                        if math.isinf(want):
                            self.assertEqual(got, want)
                        else:
                            self.assertAlmostEqual(got, want, delta=2e-3)

    def test_errors(self):
        with self.assertRaises(ph.CompareError):
            ph.dunnett([[1.0, 2.0], [3.0, 4.0]], control=2)
        with self.assertRaises(ph.CompareError):
            ph.dunnett([[1.0, 2.0], [3.0, 4.0]], alternative="both")


def dunn_numpy(groups):
    allv = np.concatenate(groups)
    r = stats.rankdata(allv)
    n = len(allv)
    _, counts = np.unique(allv, return_counts=True)
    var = n * (n + 1) / 12 - (counts ** 3 - counts).sum() / (12 * (n - 1))
    idx = np.cumsum([0] + [len(g) for g in groups])
    mr = [r[idx[i]:idx[i + 1]].mean() for i in range(len(groups))]
    out = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            z = (mr[i] - mr[j]) / math.sqrt(var * (1 / len(groups[i]) + 1 / len(groups[j])))
            out.append((z, 2 * stats.norm.sf(abs(z))))
    return out


class TestDunn(unittest.TestCase):
    def test_kruskal_matches_scipy(self):
        for groups in (make_groups((8, 10, 6), (1, 1, 1)),
                       [np.round(g) for g in make_groups((15, 12, 20, 9), (2, 2, 2, 2))]):
            h, df, p = ph.kruskal_wallis([list(g) for g in groups])
            ref = stats.kruskal(*groups)
            self.assertAlmostEqual(h, ref.statistic, delta=1e-10)
            self.assertAlmostEqual(p, ref.pvalue, delta=1e-12)
            self.assertEqual(df, len(groups) - 1)

    def test_dunn_with_ties(self):
        groups = [np.round(g) for g in make_groups((15, 12, 20, 9), (2, 2, 2, 2), 1.0)]
        ref = dunn_numpy(groups)
        for method, sm in (("bonferroni", "bonferroni"), ("sidak", "sidak"),
                           ("holm", "holm"), ("bh", "fdr_bh")):
            res = ph.dunn([list(g) for g in groups], method)
            adjusted = multipletests([p for _, p in ref], method=sm)[1]
            for c, (z, p), padj in zip(res.comparisons, ref, adjusted):
                self.assertAlmostEqual(c.statistic, z, delta=1e-10)
                self.assertAlmostEqual(c.p_raw, p, delta=1e-12)
                self.assertAlmostEqual(c.p, padj, delta=1e-12, msg=method)

    def test_p_adjust(self):
        ps = list(RNG.uniform(0, 0.2, 12)) + [1.0]
        for method, sm in (("bonferroni", "bonferroni"), ("sidak", "sidak"),
                           ("holm", "holm"), ("bh", "fdr_bh")):
            for a, b in zip(ph.p_adjust(ps, method), multipletests(ps, method=sm)[1]):
                self.assertAlmostEqual(a, b, delta=1e-12, msg=method)
        with self.assertRaises(ph.CompareError):
            ph.p_adjust(ps, "hochberg")


if __name__ == "__main__":
    unittest.main()
