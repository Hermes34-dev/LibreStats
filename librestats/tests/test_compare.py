"""Check librestats.compare (t-tests, ANOVAs) against scipy, statsmodels and numpy."""

import math
import os
import sys
import unittest

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.oneway import anova_oneway

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import compare as c  # noqa: E402

RNG = np.random.default_rng(47)
ALT = {"two-sided": "two-sided", "greater": "greater", "less": "less"}


def close(test, a, b, delta=1e-10, msg=None):
    if math.isinf(b):
        test.assertEqual(a, b, msg=msg)
    else:
        test.assertAlmostEqual(a, b, delta=delta * max(1.0, abs(b)), msg=msg)


class TestTTests(unittest.TestCase):
    def check(self, res, ref, alt, alpha, msg):
        close(self, res.statistic, ref.statistic, msg=msg)
        close(self, res.df, ref.df, msg=msg)
        close(self, res.p, ref.pvalue, delta=1e-10, msg=msg)
        ci = ref.confidence_interval(1 - alpha)
        close(self, res.ci[0], ci.low, delta=1e-8, msg=msg)
        close(self, res.ci[1], ci.high, delta=1e-8, msg=msg)

    def test_one_sample(self):
        for n in (2, 5, 30, 400):
            x = RNG.normal(10, 3, n)
            for alt in c.ALTERNATIVES:
                for alpha in (0.05, 0.01):
                    res = c.one_sample_t(list(x), 9.0, alt, alpha)
                    ref = stats.ttest_1samp(x, 9.0, alternative=alt)
                    self.check(res, ref, alt, alpha, (n, alt, alpha))
                    close(self, res.effect, (x.mean() - 9) / x.std(ddof=1))
                    close(self, res.estimate, x.mean())

    def test_welch(self):
        for n1, n2, s2 in ((5, 5, 1), (8, 30, 4), (2, 3, 0.2), (100, 60, 1.5)):
            x1, x2 = RNG.normal(5, 1, n1), RNG.normal(5.5, s2, n2)
            for alt in c.ALTERNATIVES:
                res = c.welch_t(list(x1), list(x2), alt, 0.05)
                ref = stats.ttest_ind(x1, x2, equal_var=False, alternative=alt)
                self.check(res, ref, alt, 0.05, (n1, n2, alt))
            d = (x1.mean() - x2.mean()) / math.sqrt((x1.var(ddof=1) + x2.var(ddof=1)) / 2)
            close(self, res.effect, d)
            self.assertEqual(res.n, n1 + n2)

    def test_paired(self):
        for n in (2, 6, 25, 300):
            x1 = RNG.normal(20, 4, n)
            x2 = x1 + RNG.normal(0.8, 1, n)
            for alt in c.ALTERNATIVES:
                res = c.paired_t(list(x1), list(x2), alt, 0.05)
                ref = stats.ttest_rel(x1, x2, alternative=alt)
                self.check(res, ref, alt, 0.05, (n, alt))
            d = x1 - x2
            close(self, res.effect, d.mean() / d.std(ddof=1))

    def test_paired_drops_incomplete_pairs(self):
        x1, x2 = [1.0, 2.0, None, 4.0, 6.0], [1.5, None, 3.0, 4.2, 5.0]
        ref = stats.ttest_rel([1.0, 4.0, 6.0], [1.5, 4.2, 5.0])
        res = c.paired_t(x1, x2)
        self.assertEqual(res.n, 3)
        close(self, res.statistic, ref.statistic)

    def test_by_hand(self):
        # mean 3, SD 1.5811, SE 0.7071, t = 1 / 0.7071 = 1.4142 with df 4
        res = c.one_sample_t([1.0, 2.0, 3.0, 4.0, 5.0], 2.0)
        self.assertAlmostEqual(res.statistic, math.sqrt(2), delta=1e-12)
        self.assertAlmostEqual(res.p, 2 * stats.t.sf(math.sqrt(2), 4), delta=1e-12)

    def test_errors(self):
        for call in (lambda: c.one_sample_t([1.0]), lambda: c.one_sample_t([2.0, 2.0, 2.0]),
                     lambda: c.welch_t([1.0], [1.0, 2.0]), lambda: c.welch_t([1.0, 1.0], [3.0, 3.0]),
                     lambda: c.paired_t([1.0, 2.0], [1.0]),
                     lambda: c.paired_t([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]),
                     lambda: c.one_sample_t([1.0, 2.0], alternative="both")):
            with self.assertRaises(c.CompareError):
                call()


class TestOneWay(unittest.TestCase):
    def test_matches_scipy(self):
        for sizes in ((5, 5), (8, 12, 10), (3, 20, 7, 15, 2)):
            groups = [RNG.normal(10 + i * 0.3, 1 + i * 0.2, n) for i, n in enumerate(sizes)]
            res = c.one_way_anova([list(g) for g in groups])
            ref = stats.f_oneway(*groups)
            close(self, res.statistic, ref.statistic)
            close(self, res.p, ref.pvalue)
            self.assertEqual((res.df1, res.df2), (len(sizes) - 1, sum(sizes) - len(sizes)))
            allv = np.concatenate(groups)
            ss_b = sum(len(g) * (g.mean() - allv.mean()) ** 2 for g in groups)
            ss_t = ((allv - allv.mean()) ** 2).sum()
            ms_w = (ss_t - ss_b) / res.df2
            close(self, dict(res.effects)["Eta squared (η²)"], ss_b / ss_t)
            close(self, dict(res.effects)["Omega squared (ω²)"],
                  (ss_b - res.df1 * ms_w) / (ss_t + ms_w))

    def test_welch_matches_statsmodels(self):
        for sizes, sds in (((5, 5), (1, 3)), ((8, 12, 10), (1, 1, 1)),
                           ((3, 20, 7, 15, 4), (0.5, 2, 1, 3, 1.2))):
            groups = [RNG.normal(10 + i * 0.4, s, n) for i, (n, s) in enumerate(zip(sizes, sds))]
            res = c.welch_anova([list(g) for g in groups])
            ref = anova_oneway(groups, use_var="unequal", welch_correction=True)
            close(self, res.statistic, ref.statistic)
            close(self, res.df2, ref.df[1])
            close(self, res.p, ref.pvalue)

    def test_welch_two_groups_is_welch_t(self):
        x1, x2 = list(RNG.normal(0, 1, 9)), list(RNG.normal(1, 3, 14))
        t, f = c.welch_t(x1, x2), c.welch_anova([x1, x2])
        close(self, f.statistic, t.statistic ** 2)
        close(self, f.df2, t.df)
        close(self, f.p, t.p)

    def test_errors(self):
        for call in (lambda: c.one_way_anova([[1.0, 2.0]]),
                     lambda: c.one_way_anova([[1.0], [2.0]]),
                     lambda: c.one_way_anova([[1.0, 1.0], [2.0, 2.0]]),
                     lambda: c.welch_anova([[1.0, 2.0], [3.0]]),
                     lambda: c.welch_anova([[1.0, 2.0], [3.0, 3.0]])):
            with self.assertRaises(c.CompareError):
                call()


def sphericity_numpy(x):
    """Independent check: epsilons from the eigenvalues of the double-centred
    covariance matrix, Mauchly's W from those of the contrast covariance."""
    n, k = x.shape
    s = np.cov(x, rowvar=False)
    centre = np.eye(k) - np.ones((k, k)) / k
    lam = np.linalg.eigvalsh(centre @ s @ centre)[1:]   # drop the zero eigenvalue
    gg = lam.sum() ** 2 / ((k - 1) * (lam ** 2).sum())
    hf = min(1.0, (n * (k - 1) * gg - 2) / ((k - 1) * (n - 1 - (k - 1) * gg)))
    p = k - 1
    w = lam.prod() / lam.mean() ** p
    chi2 = -(n - 1 - (2 * p * p + p + 2) / (6 * p)) * math.log(w)
    df = p * (p + 1) / 2 - 1
    return gg, hf, w, chi2, stats.chi2.sf(chi2, df)


class TestRepeatedMeasures(unittest.TestCase):
    def data(self, n, k, spread):
        subject = RNG.normal(0, 2, (n, 1))
        return subject + np.arange(k) * 0.3 + RNG.normal(0, 1, (n, k)) * np.array(spread)

    def test_matches_statsmodels(self):
        for n, k in ((4, 2), (6, 3), (12, 4), (30, 5)):
            x = self.data(n, k, np.linspace(1, 2, k))
            res = c.rm_anova(x.tolist())
            long = pd.DataFrame([(i, j, x[i, j]) for i in range(n) for j in range(k)],
                                columns=["subject", "cond", "y"])
            ref = AnovaRM(long, "y", "subject", within=["cond"]).fit().anova_table
            close(self, res.statistic, ref["F Value"].iloc[0], delta=1e-9)
            close(self, res.p, ref["Pr > F"].iloc[0], delta=1e-9)
            self.assertEqual((res.df1, res.df2), (k - 1, (k - 1) * (n - 1)))

    def test_sphericity_matches_numpy(self):
        for n, k, spread in ((6, 3, (1, 1, 3)), (12, 4, (1, 2, 1, 4)), (40, 5, (1,) * 5),
                             (8, 6, (1, 2, 3, 1, 2, 3))):
            x = self.data(n, k, spread)
            res = c.rm_anova(x.tolist())
            gg, hf, w, chi2, p = sphericity_numpy(x)
            close(self, res.sphericity["gg"], gg, delta=1e-10)
            close(self, res.sphericity["hf"], hf, delta=1e-10)
            close(self, res.sphericity["W"], w, delta=1e-9)
            close(self, res.sphericity["chi2"], chi2, delta=1e-9)
            close(self, res.sphericity["p"], p, delta=1e-9)
            names = [r[0] for r in res.corrections]
            self.assertEqual(names, ["Sphericity assumed", "Greenhouse-Geisser", "Huynh-Feldt",
                                     "Lower bound"])
            name, eps, df1, df2, pc = res.corrections[1]
            close(self, pc, stats.f.sf(res.statistic, gg * (k - 1), gg * (k - 1) * (n - 1)),
                  delta=1e-9)

    def test_two_conditions_is_paired_t(self):
        x = self.data(10, 2, (1, 2))
        res = c.rm_anova(x.tolist())
        t = c.paired_t(list(x[:, 0]), list(x[:, 1]))
        close(self, res.statistic, t.statistic ** 2)
        close(self, res.p, t.p)
        self.assertEqual(res.sphericity["gg"], 1.0)
        self.assertNotIn("W", res.sphericity)

    def test_effects_and_incomplete_rows(self):
        x = self.data(9, 3, (1, 1, 1)).tolist()
        res = c.rm_anova(x + [[1.0, None, 2.0]])
        self.assertEqual(res.n, 9)
        ss = {r[0]: r[1] for r in res.table}
        close(self, dict(res.effects)["Partial eta squared (η²p)"],
              ss["Conditions"] / (ss["Conditions"] + ss["Error"]))
        close(self, ss["Total"], ss["Conditions"] + ss["Subjects"] + ss["Error"])

    def test_mauchly_needs_enough_subjects(self):
        res = c.rm_anova(self.data(3, 4, (1, 1, 1, 1)).tolist())
        self.assertNotIn("W", res.sphericity)
        self.assertTrue(res.notes)

    def test_errors(self):
        for rows in ([[1.0, 2.0]], [[1.0], [2.0]], [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]):
            with self.assertRaises(c.CompareError, msg=rows):
                c.rm_anova(rows)


if __name__ == "__main__":
    unittest.main()
