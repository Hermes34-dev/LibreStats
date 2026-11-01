"""Check librestats.screening against Montgomery's worked examples and statsmodels OLS."""

import math
import os
import sys
import unittest
from itertools import combinations

import numpy as np
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import doe  # noqa: E402
from librestats import screening as sc  # noqa: E402

RNG = np.random.default_rng(21)

# Montgomery, Design and Analysis of Experiments, Example 6.2: filtration rate, 2^4,
# standard order
FILTRATION = [45, 71, 48, 65, 68, 60, 80, 65, 43, 100, 45, 104, 75, 86, 70, 96]
FILTRATION_EFFECTS = {"A": 21.625, "B": 3.125, "C": 9.875, "D": 14.625, "AB": 0.125,
                      "AC": -18.125, "AD": 16.625, "BC": 2.375, "BD": -0.375, "CD": -1.125,
                      "ABC": 1.875, "ABD": 4.125, "ACD": -1.625, "BCD": -2.625, "ABCD": 1.375}


def ols_reference(coded, ptype, blocks, y, masks):
    """statsmodels fit of the same model: intercept, block dummies, centre term, effects."""
    obs = [i for i, v in enumerate(y) if v is not None]
    levels = sorted({blocks[i] for i in obs})
    centre = any(ptype[i] == 0 for i in obs)
    X = []
    for i in obs:
        row = [1.0] + [float(blocks[i] == b) for b in levels[1:]]
        if centre:
            row.append(float(ptype[i] == 0))
        row += [float(np.prod([coded[i][j] for j in range(len(coded[i])) if m >> j & 1]))
                for m in masks]
        X.append(row)
    return sm.OLS(np.array([y[i] for i in obs]), np.array(X)).fit(), len(X[0]) - len(masks)


class TestMontgomery(unittest.TestCase):
    def test_filtration_effects_and_lenth(self):
        a = sc.analyse(doe.standard_order(4), [1] * 16, [1] * 16, list(map(float, FILTRATION)))
        got = {t.name: e for t, e in a.contrasts}
        self.assertEqual(set(got), set(FILTRATION_EFFECTS))
        for name, effect in FILTRATION_EFFECTS.items():
            self.assertAlmostEqual(got[name], effect, places=10)
        # Montgomery (section 6.5): PSE = 2.625, ME = 6.75, SME = 13.70
        self.assertAlmostEqual(a.lenth.pse, 2.625)
        self.assertAlmostEqual(a.lenth.me, 6.75, places=2)
        self.assertAlmostEqual(a.lenth.sme, 13.70, places=2)
        self.assertEqual([t.name for t in a.terms],
                         ["A", "B", "C", "D", "AB", "AC", "AD", "BC", "BD", "CD"])
        self.assertEqual((a.df, a.error), (5, "residual"))

    def test_saturated_model_uses_lenth(self):
        # 2^(3-1): A = BC, B = AC, C = AB; three terms, four runs, no residual df
        x = doe.fractional_factorial(3, 4).matrix()
        a = sc.analyse(x, [1] * 4, [1] * 4, [10.0, 14.0, 13.0, 30.0])
        self.assertEqual((a.df, a.error), (0, "lenth"))
        self.assertEqual([(t.name, t.alias_text()) for t in a.terms],
                         [("A", "BC"), ("B", "AC"), ("C", "AB")])
        for t in a.terms:
            self.assertAlmostEqual(t.t, t.effect / a.lenth.pse)
            self.assertAlmostEqual(t.p, 2 * stats.t.sf(abs(t.t), a.lenth.m / 3))

    def test_curvature_example_6_7(self):
        # 2^2 with 5 centre points: yF = 40.425, yC = 40.46, SS = 0.0027, MSE = 0.0430
        x = doe.standard_order(2) + [[0, 0]] * 5
        y = [39.3, 40.9, 40.0, 41.5, 40.3, 40.5, 40.7, 40.2, 40.6]
        a = sc.analyse(x, [1] * 4 + [0] * 5, [1] * 9, y)
        c = a.curvature
        self.assertAlmostEqual(c.mean_factorial, 40.425)
        self.assertAlmostEqual(c.mean_centre, 40.46)
        self.assertAlmostEqual(c.ss, 4 * 5 * 0.035 ** 2 / 9)
        self.assertAlmostEqual(a.pure_ss / a.pure_df, 0.043)
        self.assertAlmostEqual(c.f, c.ss / 0.043)
        self.assertAlmostEqual(c.p, stats.f.sf(c.f, 1, 4))
        self.assertEqual([round(t.effect, 10) for t in a.terms], [1.55, 0.65, -0.05])


class TestAgainstStatsmodels(unittest.TestCase):
    def check(self, coded, ptype, blocks, y, model="2fi"):
        a = sc.analyse(coded, ptype, blocks, y, model)
        fit, first = ols_reference(coded, ptype, blocks, y, [t.mask for t in a.terms])
        self.assertEqual(a.df, int(fit.df_resid))
        np.testing.assert_allclose([t.coef for t in a.terms], fit.params[first:], atol=1e-9)
        np.testing.assert_allclose([t.se / 2 for t in a.terms], fit.bse[first:], rtol=1e-8)
        np.testing.assert_allclose([t.p for t in a.terms], fit.pvalues[first:], rtol=1e-6,
                                   atol=1e-12)
        self.assertAlmostEqual(a.r2, fit.rsquared, places=10)
        self.assertAlmostEqual(a.r2_adj, fit.rsquared_adj, places=10)
        return a

    def test_replicated_blocked_with_centre_points_and_missing(self):
        d = doe.build_design([doe.Factor(c, "", -1.0, 1.0) for c in "PQRS"], "full",
                             centre_points=3, replicates=2, blocks=2, seed=4)
        coded = [r.coded for r in d.runs]
        ptype = [0 if r.centre else 1 for r in d.runs]
        blocks = [r.block for r in d.runs]
        y = [float(5 + 3 * c[0] - 2 * c[1] * c[2] + r.block + RNG.normal()) for c, r in
             zip(coded, d.runs)]
        y[3] = y[17] = None
        a = self.check(coded, ptype, blocks, y)
        self.assertEqual(a.missing, 2)
        self.assertIsNotNone(a.curvature)
        self.assertEqual(a.contrasts, [])          # incomplete data: no Lenth

    def test_fraction_main_effects_only(self):
        x = doe.fractional_factorial(7, 16).matrix()
        y = list(RNG.normal(size=16))
        a = self.check(x, [1] * 16, [1] * 16, y, model="main")
        self.assertEqual([t.name for t in a.terms], list("ABCDEFG"))


class TestAliasing(unittest.TestCase):
    def test_chains_follow_the_defining_relation(self):
        f = doe.fractional_factorial(7, 16)
        a = sc.analyse(f.matrix(), [1] * 16, [1] * 16, list(RNG.normal(size=16)))
        self.assertEqual(len(a.contrasts), 15)
        for t, _ in a.contrasts:
            members = {t.mask} | {m for _, m in t.aliases}
            self.assertEqual(members, f.aliases(t.mask))
        # 2fi model: 7 main effects + the 7 two-factor chains (21 2fis in 7 chains)
        self.assertEqual(len(a.terms), 14)
        self.assertEqual(a.terms[7].name, "AB")
        self.assertEqual(a.terms[7].alias_text(max_order=2), "CE + DF (+ 5 higher-order)")

    def test_negative_aliases(self):
        x = [[a, b, -a * b] for b in (-1, 1) for a in (-1, 1)]       # C = -AB
        a = sc.analyse(x, [1] * 4, [1] * 4, [1.0, 2.0, 4.0, 3.0], model="main")
        self.assertEqual(a.terms[2].alias_text(), "- AB")

    def test_blocks_confounded_effect_is_excluded(self):
        d = doe.build_design([doe.Factor(c, "", -1.0, 1.0) for c in "ABCD"], "full",
                             blocks=2, randomize=False)
        a = sc.analyse([r.coded for r in d.runs], [1] * 16, [r.block for r in d.runs],
                       list(map(float, FILTRATION)))
        self.assertEqual([t.name for t in a.block_effects], ["ABCD"])
        self.assertEqual(len(a.contrasts), 14)

    def test_plackett_burman_partial_aliasing(self):
        x = doe.plackett_burman(11, 12)
        a = sc.analyse(x, [1] * 12, [1] * 12, list(RNG.normal(size=12)))
        self.assertEqual([t.name for t in a.terms], list("ABCDEFGHJKL"))
        self.assertTrue(any("partially aliased" in n for n in a.notes))
        self.assertEqual((a.df, a.error), (0, "lenth"))


class TestLenth(unittest.TestCase):
    def test_against_direct_formula(self):
        for m in (7, 15, 31):
            c = RNG.normal(size=m) * np.where(RNG.random(m) < 0.2, 6, 1)
            a = np.abs(c)
            s0 = 1.5 * np.median(a)
            pse = 1.5 * np.median(a[a < 2.5 * s0])
            d = m / 3
            gamma = (1 + 0.95 ** (1 / m)) / 2
            L = sc.lenth(list(c))
            self.assertAlmostEqual(L.pse, pse, places=12)
            self.assertAlmostEqual(L.me, stats.t.ppf(0.975, d) * pse, places=6)
            self.assertAlmostEqual(L.sme, stats.t.ppf(gamma, d) * pse, places=6)

    def test_half_normal_quantiles(self):
        q = sc.half_normal_quantiles(15)
        ref = stats.halfnorm.ppf((np.arange(1, 16) - 0.5) / 15)
        np.testing.assert_allclose(q, ref, atol=1e-9)


class TestErrors(unittest.TestCase):
    def test_refusals(self):
        x = doe.standard_order(2)
        with self.assertRaisesRegex(sc.ScreeningError, "response surface"):
            sc.analyse(x + [[1.4, 0]], [1, 1, 1, 1, -1], [1] * 5, [1.0] * 5)
        with self.assertRaisesRegex(sc.ScreeningError, "at least 3"):
            sc.analyse(x, [1] * 4, [1] * 4, [1.0, None, None, 2.0])
        with self.assertRaisesRegex(sc.ScreeningError, "-1 or \\+1"):
            sc.analyse([[0.5, 1]] + x[1:], [1] * 4, [1] * 4, [1.0] * 4)
        a = sc.analyse(x + [[0, 0]], [1] * 4 + [0], [1] * 5, [1.0, 2.0, 3.0, 5.0, 2.0])
        self.assertIsNone(a.curvature.f)
        self.assertTrue(any("at least 2 centre points" in n for n in a.notes))


if __name__ == "__main__":
    unittest.main()
