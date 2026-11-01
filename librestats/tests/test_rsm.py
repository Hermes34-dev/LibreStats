"""Check librestats.rsm against Montgomery's Example 11.2, statsmodels and numpy."""

import os
import sys
import unittest
from itertools import combinations, product

import numpy as np
import statsmodels.api as sm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import doe, rsm  # noqa: E402

RNG = np.random.default_rng(8)
A = 2 ** 0.5

# Montgomery, Design and Analysis of Experiments, Example 11.2: CCD, time (80-90 min) and
# temperature (170-180 F), yield
CCD_X = [[-1, -1], [-1, 1], [1, -1], [1, 1]] + [[0, 0]] * 5 + [[A, 0], [-A, 0], [0, A], [0, -A]]
CCD_Y = [76.5, 77.0, 78.0, 79.5, 79.9, 80.3, 80.0, 79.7, 79.8, 78.4, 75.6, 78.5, 77.0]


def quad(x):
    k = len(x)
    return [1.0] + list(x) + [v * v for v in x] + [x[i] * x[j] for i, j in combinations(range(k), 2)]


class TestMontgomery(unittest.TestCase):
    def setUp(self):
        self.fit = rsm.fit_model(CCD_X, [1] * 13, CCD_Y)

    def test_coefficients(self):
        # y = 79.94 + 0.995 x1 + 0.515 x2 - 1.376 x1^2 - 1.001 x2^2 + 0.250 x1 x2
        ref = {"Constant": 79.94, "A": 0.995, "B": 0.515, "A^2": -1.376, "B^2": -1.001,
               "AB": 0.25}
        for t in self.fit.terms:
            self.assertAlmostEqual(t.coef, ref[t.name], delta=0.0006)

    def test_anova_lack_of_fit(self):
        rows = {r.source.strip(): r for r in self.fit.anova}
        self.assertAlmostEqual(rows["Residual error"].ss, 0.4953, places=4)
        self.assertAlmostEqual(rows["Lack of fit"].ss, 0.2833, places=4)
        self.assertAlmostEqual(rows["Pure error"].ss, 0.2120, places=4)
        self.assertEqual((rows["Lack of fit"].df, rows["Pure error"].df), (3, 4))

    def test_stationary_point_and_canonical_analysis(self):
        st = self.fit.stationary
        np.testing.assert_allclose(st.coded, [0.389, 0.306], atol=0.0006)
        self.assertAlmostEqual(st.response, 80.21, places=2)
        np.testing.assert_allclose(st.eigenvalues, [-1.4141, -0.9634], atol=0.0001)
        self.assertEqual((st.nature, st.inside), ("maximum", True))
        c, h, _ = rsm.coding(CCD_X, [[85 + 5 * a, 175 + 5 * b] for a, b in CCD_X])
        self.assertAlmostEqual(c[0] + h[0] * st.coded[0], 86.95, places=2)    # minutes
        self.assertAlmostEqual(c[1] + h[1] * st.coded[1], 176.53, places=2)   # degrees F

    def test_equation_in_actual_units(self):
        # -1430.52 + 7.81 t + 13.27 T - 0.055 t^2 - 0.040 T^2 + 0.010 t T
        eq = dict(rsm.actual_equation(self.fit, [85, 175], [5, 5]))
        for name, value, tol in (("Constant", -1430.52, 0.01), ("A", 7.81, 0.005),
                                 ("B", 13.27, 0.005), ("A^2", -0.055, 0.0006),
                                 ("B^2", -0.040, 0.0006), ("AB", 0.010, 0.0001)):
            self.assertAlmostEqual(eq[name], value, delta=tol)

    def test_maximum_within_the_region_is_the_stationary_point(self):
        np.testing.assert_allclose(self.fit.maximum[0], self.fit.stationary.coded, atol=1e-9)
        self.assertEqual(self.fit.bounds, [(-A, A), (-A, A)])


class TestAgainstStatsmodels(unittest.TestCase):
    def test_box_behnken_with_blocks_and_missing(self):
        d = doe.build_rsm([doe.Factor(f"X{i}", "", 0.0, 1.0) for i in range(4)], "bbd",
                          blocks=3, seed=2)
        x = [r.coded for r in d.runs]
        blocks = [r.block for r in d.runs]
        y = [float(10 + v[0] - 2 * v[1] ** 2 + v[2] * v[3] + b + RNG.normal(0, 0.3))
             for v, b in zip(x, blocks)]
        y[5] = None
        fit = rsm.fit_model(x, blocks, y)
        obs = [i for i, v in enumerate(y) if v is not None]
        Z = [[1.0 if blocks[i] == b else -1.0 if blocks[i] == 3 else 0.0 for b in (1, 2)]
             for i in obs]
        X = np.array([[1.0] + z + quad(x[i])[1:] for z, i in zip(Z, obs)])
        ref = sm.OLS(np.array([y[i] for i in obs]), X).fit()
        np.testing.assert_allclose([t.coef for t in fit.terms], ref.params, atol=1e-9)
        np.testing.assert_allclose([t.se for t in fit.terms], ref.bse, rtol=1e-8)
        np.testing.assert_allclose([t.p for t in fit.terms], ref.pvalues, rtol=1e-6,
                                   atol=1e-14)
        self.assertAlmostEqual(fit.r2_adj, ref.rsquared_adj, places=10)
        h = ref.get_influence().hat_matrix_diag
        press = np.sum((ref.resid / (1 - h)) ** 2)
        self.assertAlmostEqual(fit.r2_pred, 1 - press / ref.centered_tss, places=10)
        # group F test = statsmodels' F test of the same restriction
        square = [j for j, t in enumerate(fit.terms) if t.kind == "square"]
        R = np.zeros((len(square), len(fit.terms)))
        for r, j in enumerate(square):
            R[r, j] = 1
        row = next(r for r in fit.anova if r.source == "Square")
        self.assertAlmostEqual(row.f, float(ref.f_test(R).fvalue), places=8)

    def test_linear_and_interaction_models(self):
        x = doe.standard_order(3) + [[0, 0, 0]] * 4
        y = list(RNG.normal(size=12))
        fit = rsm.fit_model(x, [1] * 12, y, "interaction")
        self.assertEqual([t.name for t in fit.terms],
                         ["Constant", "A", "B", "C", "AB", "AC", "BC"])
        self.assertIsNone(fit.stationary)
        with self.assertRaisesRegex(rsm.RSMError, "cannot separate the square terms"):
            rsm.fit_model(x, [1] * 12, y, "quadratic")


class TestAlgebra(unittest.TestCase):
    def test_eigen_against_numpy(self):
        for k in (2, 3, 5, 7):
            m = RNG.normal(size=(k, k))
            m = (m + m.T) / 2
            values, vectors = rsm.eigen_symmetric(m.tolist())
            np.testing.assert_allclose(values, np.linalg.eigvalsh(m), atol=1e-10)
            for v, vec in zip(values, vectors):
                np.testing.assert_allclose(m @ vec, v * np.array(vec), atol=1e-9)

    def test_box_optimum_against_brute_force(self):
        for seed in range(4):
            rng = np.random.default_rng(seed)
            k = 3
            x = [list(v) for v in product((-1, 0, 1), repeat=k)]
            y = [float(rng.normal()) + sum(rng.normal() * a * b for a, b in combinations(p, 2))
                 for p in x]
            fit = rsm.fit_model(x, [1] * len(x), y)
            grid = np.linspace(-1, 1, 41)
            values = [(fit.predict(list(p)), p) for p in product(grid, repeat=k)]
            best = max(values)[0]
            worst = min(values)[0]
            self.assertGreaterEqual(fit.maximum[1], best - 1e-12)
            self.assertLessEqual(fit.minimum[1], worst + 1e-12)
            self.assertLess(fit.maximum[1] - best, 0.01)    # grid is close to the exact optimum
            self.assertTrue(all(-1 - 1e-12 <= v <= 1 + 1e-12 for v in fit.maximum[0]))

    def test_saddle_point(self):
        x = [list(v) for v in product((-1, 0, 1), repeat=2)]
        y = [a * a - b * b + 0.1 * a for a, b in x]
        st = rsm.fit_model(x, [1] * 9, y).stationary
        self.assertEqual(st.nature, "saddle point")

    def test_actual_equation_reproduces_predictions(self):
        x = [list(v) for v in product((-1, 0, 1), repeat=3)]
        y = [float(v) for v in RNG.normal(size=27)]
        fit = rsm.fit_model(x, [1] * 27, y)
        c, h = [60.0, 5.0, 0.3], [20.0, 2.5, 0.1]
        eq = dict(rsm.actual_equation(fit, c, h))
        for p in RNG.uniform(-1.5, 1.5, (5, 3)):
            X = [ci + hi * pi for ci, hi, pi in zip(c, h, p)]
            value = eq["Constant"] + sum(eq[n] * X[i] for i, n in enumerate("ABC"))
            value += sum(eq[n + "^2"] * X[i] ** 2 for i, n in enumerate("ABC"))
            value += sum(eq[a + b] * X[i] * X[j] for (i, a), (j, b) in
                         combinations(enumerate("ABC"), 2))
            self.assertAlmostEqual(value, fit.predict(list(p)), places=8)


class TestCoding(unittest.TestCase):
    def test_exact_from_rounded_design_sheets(self):
        factors = [doe.Factor("T", "", 40.0, 80.0), doe.Factor("p", "", 0.15, 0.35),
                   doe.Factor("s", "", 200.0, 600.0)]
        for kind, kw in (("ccd", {}), ("ccd", {"inscribed": True}), ("ccd", {"alpha": "face"}),
                         ("bbd", {})):
            d = doe.build_rsm(factors, kind, randomize=False, **kw)
            coded = [r.coded for r in d.runs]
            actual = [[f.setting(v, d.scale) for f, v in zip(factors, r.coded)] for r in d.runs]
            c, h, _ = rsm.coding(coded, actual)
            for f, ci, hi in zip(factors, c, h):
                self.assertAlmostEqual(ci, (f.low + f.high) / 2, places=12)
                self.assertAlmostEqual(hi * d.scale, (f.high - f.low) / 2, places=12)


if __name__ == "__main__":
    unittest.main()
