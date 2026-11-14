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

    def test_r2_ceiling(self):
        # Montgomery: pure error SS 0.2120, total SS 28.7431; the ceiling bounds R^2
        self.assertAlmostEqual(self.fit.r2_max, 1 - 0.2120 / 28.743077, places=6)
        self.assertLess(self.fit.r2, self.fit.r2_max)
        # independent: 1 - within-group SS of repeated points / total SS (numpy)
        y = np.array(CCD_Y)
        centre = y[4:9]
        self.assertAlmostEqual(self.fit.r2_max, 1 - ((centre - centre.mean()) ** 2).sum() /
                               ((y - y.mean()) ** 2).sum(), places=12)
        # no repeated runs: no ceiling
        fit = rsm.fit_model(CCD_X[:5] + CCD_X[9:], [1] * 9, CCD_Y[:5] + CCD_Y[9:])
        self.assertIsNone(fit.r2_max)
        # blocks: pure error only from repeats within the same block
        fit = rsm.fit_model(CCD_X, [1] * 7 + [2] * 6, CCD_Y)
        g1, g2 = y[4:7], y[7:9]
        pe = ((g1 - g1.mean()) ** 2).sum() + ((g2 - g2.mean()) ** 2).sum()
        self.assertAlmostEqual(fit.r2_max, 1 - pe / ((y - y.mean()) ** 2).sum(), places=12)

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


class TestCCCDAnalysis(unittest.TestCase):
    """A CCCD sheet (Low/High = +-alpha) analysed as any response surface design."""

    def test_coding_and_fit(self):
        factors = [doe.Factor("T", "C", 30.0, 90.0), doe.Factor("t", "min", 2.5, 22.5),
                   doe.Factor("S", "%", 0.0, 100.0)]
        d = doe.build_rsm(factors, "cccd", seed=11)
        coded = [r.coded for r in d.runs]
        actual = [[f.setting(c, d.scale) for f, c in zip(factors, r.coded)] for r in d.runs]
        centres, halves, worst = rsm.coding(coded, actual)
        for f, c, h in zip(factors, centres, halves):
            self.assertAlmostEqual(c, (f.low + f.high) / 2, places=12)
            self.assertAlmostEqual(h, (f.high - f.low) / (2 * d.alpha), places=12)
        beta = RNG.normal(size=10)
        y = [float(np.dot(beta, quad(x))) + RNG.normal(scale=0.1) for x in coded]
        fit = rsm.fit_model(coded, [r.block for r in d.runs], y)
        ref = sm.OLS(y, np.array([quad(x) for x in coded])).fit()
        self.assertTrue(np.allclose([t.coef for t in fit.terms], ref.params, atol=1e-9))
        rows = {r.source.strip(): r for r in fit.anova}
        self.assertEqual((rows["Residual error"].df, rows["Lack of fit"].df,
                          rows["Pure error"].df), (18, 13, 5))


def cubic_row(x):
    """Full cubic model row (numpy reference): all monomials of degree <= 3."""
    from itertools import combinations_with_replacement as cwr
    return [float(np.prod([x[i] for i in m])) for d in range(4) for m in cwr(range(len(x)), d)]


def standard_designs():
    f = lambda k: [doe.Factor(f"X{i}", "", 10.0, 30.0) for i in range(k)]  # noqa: E731
    return [("CCD 2", doe.build_rsm(f(2), "ccd", randomize=False)),
            ("CCD 3", doe.build_rsm(f(3), "ccd", randomize=False)),
            ("CCD 4 face", doe.build_rsm(f(4), "ccd", alpha="face", randomize=False)),
            ("BBD 3", doe.build_rsm(f(3), "bbd", randomize=False)),
            ("BBD 4 blocks", doe.build_rsm(f(4), "bbd", blocks=3, randomize=False)),
            ("CCCD 3", doe.build_rsm(f(3), "cccd", randomize=False)),
            ("CCD 3 blocks", doe.build_rsm(f(3), "ccd", blocks=3, randomize=False)),
            ("2^3 + centre", doe.build_design(f(3), "full", centre_points=3,
                                              randomize=False))]


class TestReducedCubic(unittest.TestCase):
    def fit(self, d, seed=0):
        x = [r.coded for r in d.runs]
        rng = np.random.default_rng(seed)
        y = [float(5 + v[0] - v[1] ** 2 + 0.4 * v[0] * v[1] * v[-1] + 0.3 * v[0] ** 3 +
                   rng.normal(scale=0.2)) for v in x]
        return x, y, rsm.fit_model(x, [r.block for r in d.runs], y, "cubic")

    def test_rank_aliases_and_fit_against_numpy(self):
        for name, d in standard_designs():
            x, y, fit = self.fit(d)
            blocks = [r.block for r in d.runs]
            levels = sorted(set(blocks))
            B = np.array([[1.0 if b == lv else -1.0 if b == levels[-1] else 0.0
                           for lv in levels[:-1]] for b in blocks]).reshape(len(x), -1)
            full = np.hstack([np.array([cubic_row(v) for v in x]), B])
            rank = np.linalg.matrix_rank(full)
            self.assertEqual(len(fit.terms), rank, name)
            self.assertEqual(fit.full_terms, full.shape[1], name)
            self.assertEqual(len(fit.aliases), full.shape[1] - rank, name)
            # same column space as the full cubic: same fitted values and residual SS
            beta, *_ = np.linalg.lstsq(full, y, rcond=None)
            sse_ref = float(((np.array(y) - full @ beta) ** 2).sum())
            self.assertAlmostEqual(fit.mse * fit.df, sse_ref, places=8, msg=name)
            np.testing.assert_allclose([g.fitted for g in fit.diagnostics], full @ beta,
                                       atol=1e-8, err_msg=name)
            # every alias is exact: the removed column equals the stated combination
            col = {t.name: [rsm._value(t, v) for v in x] for t in rsm.model_terms(len(x[0]),
                                                                                  "cubic")}
            for a in fit.aliases:
                combo = sum(c * np.array(col[n]) for c, n in a.combination) \
                    if a.combination else np.zeros(len(x))
                np.testing.assert_allclose(col[a.term], combo, atol=1e-9, err_msg=a.text())

    def test_hierarchy_kept_in_response_surface_designs(self):
        for name, d in standard_designs()[:-1]:
            _, _, fit = self.fit(d)
            k = len(d.factors)
            self.assertEqual(sum(t.kind != "cubic" for t in fit.terms),
                             doe.model_terms(k) + len(set(r.block for r in d.runs)) - 1, name)

    def test_alias_texts(self):
        designs = dict(standard_designs())
        texts = {name: [a.text() for a in self.fit(designs[name])[2].aliases]
                 for name in designs}
        self.assertIn("A^3 is aliased with A and AB^2: A^3 = 2 A - AB^2", texts["CCD 2"])
        self.assertIn("AC^2 is aliased with AB^2", texts["CCD 3"])
        self.assertIn("A^3 is aliased with A", texts["CCD 4 face"])
        self.assertIn("ABC is zero in every run (it cannot be estimated)", texts["BBD 3"])
        self.assertIn("B^2 is aliased with A^2", texts["2^3 + centre"])
        self.assertIn("A^2B is aliased with B", texts["2^3 + centre"])
        self.assertEqual(len(texts["CCCD 3"]), 3)                  # rank 17 of 20

    def test_full_cubic_estimable(self):
        d = doe.build_rsm([doe.Factor("X0", "", 10.0, 30.0), doe.Factor("X1", "", 1.0, 2.0)],
                          "cccd", randomize=False)                  # k = 2: rank 10 of 10
        _, _, fit = self.fit(d)
        self.assertEqual((len(fit.terms), fit.aliases), (10, []))

    def test_vif_against_statsmodels(self):
        from statsmodels.stats.outliers_influence import variance_inflation_factor
        x, y, fit = self.fit(dict(standard_designs())["CCCD 3"])
        X = np.array([[rsm._value(t, v) for t in fit.terms] for v in x])
        for j, t in enumerate(fit.terms):
            if t.kind != "constant":
                self.assertAlmostEqual(t.vif, variance_inflation_factor(X, j), places=8)

    def test_numeric_optimum_against_dense_grid(self):
        for name in ("CCD 2", "CCD 3", "CCCD 3"):
            x, y, fit = self.fit(dict(standard_designs())[name], seed=4)
            k = fit.k
            m = 201 if k == 2 else 41
            axes = [np.linspace(lo, hi, m) for lo, hi in fit.bounds]
            grid = [list(p) for p in product(*axes)]
            values = [fit.predict(p) for p in grid]
            hi_x, hi_y = fit.maximum
            lo_x, lo_y = fit.minimum
            self.assertGreaterEqual(hi_y, max(values) - 1e-9, name)
            self.assertLessEqual(lo_y, min(values) + 1e-9, name)
            self.assertAlmostEqual(fit.predict(hi_x), hi_y, places=12)
            self.assertTrue(all(lo - 1e-12 <= v <= hi + 1e-12
                                for v, (lo, hi) in zip(hi_x, fit.bounds)))
            self.assertIsNone(fit.stationary)

    def test_actual_equation_reproduces_cubic_predictions(self):
        x, y, fit = self.fit(dict(standard_designs())["CCCD 3"])
        c, h = [20.0, 5.0, 80.0], [4.0, 0.5, 12.0]
        eq = rsm.actual_equation(fit, c, h)
        for v in RNG.uniform(-1.6, 1.6, size=(6, 3)):
            X = [ci + hi * vi for ci, hi, vi in zip(c, h, v)]
            pred = sum(coef * np.prod([X[i] for i in rsm.name_monomial(n)]) for n, coef in eq)
            self.assertAlmostEqual(pred, fit.predict(list(v)), places=8)


class TestMapPairs(unittest.TestCase):
    def fit(self, k, seed=2):
        d = doe.build_rsm([doe.Factor(f"X{i}", "", 0.0, 1.0) for i in range(k)], "ccd",
                          randomize=False)
        x = [r.coded for r in d.runs]
        rng = np.random.default_rng(seed)
        beta = rng.normal(size=len(quad(x[0])))
        y = [float(np.dot(beta, quad(v)) + rng.normal(scale=0.1)) for v in x]
        return rsm.fit_model(x, [1] * len(x), y)

    def test_all_pairs_up_to_three(self):
        self.assertEqual(rsm.map_pairs(self.fit(2)), ([(0, 1)], [(0, 1)]))
        pairs, every = rsm.map_pairs(self.fit(3))
        self.assertEqual(pairs, [(0, 1), (0, 2), (1, 2)])
        self.assertEqual(pairs, every)

    def test_capped_at_three_by_predicted_range(self):
        for k in (4, 5):
            fit = self.fit(k)
            pairs, every = rsm.map_pairs(fit)
            self.assertEqual(len(every), k * (k - 1) // 2)
            self.assertEqual(len(pairs), 3)
            spread = {}
            for i, j in every:                  # independent: numpy grid, others at 0
                (x0, x1), (y0, y1) = fit.bounds[i], fit.bounds[j]
                vals = []
                for a in np.linspace(x0, x1, 21):
                    for b in np.linspace(y0, y1, 21):
                        x = [0.0] * k
                        x[i], x[j] = a, b
                        vals.append(fit.predict(x))
                spread[(i, j)] = max(vals) - min(vals)
            top = sorted(sorted(every, key=lambda p: -spread[p])[:3])
            self.assertEqual(pairs, top)
