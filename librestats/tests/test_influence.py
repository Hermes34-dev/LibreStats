"""Check librestats.influence (Cook's distance) against Cook (1977) Tables 1 and 2,
statsmodels' OLSInfluence and refitting with each observation deleted (eq. 4), and
the diagnostics of the response surface and screening analyses."""

import os
import sys
import unittest

import numpy as np
import statsmodels.api as sm
from scipy.stats import f as f_dist
from statsmodels.stats.outliers_influence import OLSInfluence

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import doe, influence, rsm, screening  # noqa: E402

# Hald cement data (Draper & Smith, 1966, p. 165), Cook's Example 2
HALD_X = np.array([[7, 26, 6, 60], [1, 29, 15, 52], [11, 56, 8, 20], [11, 31, 8, 47],
                   [7, 52, 6, 33], [11, 55, 9, 22], [3, 71, 17, 6], [1, 31, 22, 44],
                   [2, 54, 18, 22], [21, 47, 4, 26], [1, 40, 23, 34], [11, 66, 9, 12],
                   [10, 68, 8, 12]], float)
HALD_Y = np.array([78.5, 74.3, 104.3, 87.6, 95.9, 109.2, 102.7, 72.5, 93.1, 115.9, 83.8,
                   113.3, 109.4])
# Cook (1977) Table 2: |t_i|, V(yhat)/V(R), D_i (None: printed as "< 2e-3")
HALD_TABLE = [(0.003, 1.22, None), (0.76, 0.50, 0.06), (1.05, 1.36, 0.30),
              (0.84, 0.42, 0.06), (0.13, 0.56, None), (1.71, 0.14, 0.08),
              (0.74, 0.58, 0.06), (1.69, 0.69, 0.31), (0.67, 0.42, 0.04),
              (0.21, 2.34, 0.02), (1.07, 0.74, 0.17), (0.46, 0.36, 0.02), (1.12, 0.44, 0.11)]
# Cook (1977) Table 1, Longley data: |t_i|, V(yhat)/V(R), D_i (None: "< 5e-3")
LONGLEY_TABLE = [(1.15, 0.74, 0.14), (0.48, 1.30, 0.04), (0.19, 0.57, None),
                 (1.70, 0.59, 0.24), (1.64, 1.60, 0.61), (1.03, 0.59, 0.09),
                 (0.75, 0.97, 0.08), (0.06, 1.02, None), (0.07, 0.84, None),
                 (1.83, 0.49, 0.23), (0.07, 0.56, None), (0.18, 0.93, None),
                 (0.64, 0.60, 0.04), (0.32, 0.30, None), (1.42, 0.59, 0.17),
                 (1.21, 2.21, 0.47)]


def ours(X, y):
    X = sm.add_constant(np.asarray(X, float))
    inv = np.linalg.inv(X.T @ X)
    beta = inv @ X.T @ y
    sse = float(((y - X @ beta) ** 2).sum())
    return X, influence.diagnostics(X.tolist(), list(map(float, y)), beta.tolist(),
                                    inv.tolist(), sse)


class TestCookPaper(unittest.TestCase):
    def check_table(self, X, y, table, tol_t, tol_ratio, tol_d, skip=()):
        _, diags = ours(X, y)
        for i, (d, (t, ratio, cooks)) in enumerate(zip(diags, table)):
            self.assertAlmostEqual(abs(d.studentized), t, delta=tol_t, msg=i)
            self.assertAlmostEqual(d.ratio, ratio, delta=tol_ratio, msg=i)
            if i in skip:
                continue
            if cooks is None:
                self.assertLess(d.cooks, 0.005, i)
            else:
                self.assertAlmostEqual(d.cooks, cooks, delta=tol_d, msg=i)
        return diags

    def test_hald_table_2(self):
        # observation 8 is printed as D = 0.31, but the same row's |t| = 1.69 and
        # V(yhat)/V(R) = 0.69 give 1.69^2 x 0.69 / 5 = 0.394 by Cook's eq. (8): a misprint
        diags = self.check_table(HALD_X, HALD_Y, HALD_TABLE, 0.006, 0.006, 0.006, skip={7})
        self.assertAlmostEqual(1.69 ** 2 * 0.69 / 5, diags[7].cooks, delta=0.005)
        self.assertEqual(max(range(13), key=lambda i: diags[i].cooks), 7)  # largest: obs 8
        self.assertFalse(any(d.influential for d in diags))   # "data seem well behaved"

    def test_longley_table_1(self):
        # Longley's data are badly conditioned; the 1977 printout differs by up to 0.013
        # in |t| (1948: 0.468 here and in statsmodels, 0.48 printed)
        lo = sm.datasets.longley.load_pandas().data
        diags = self.check_table(lo.iloc[:, 1:].values, lo.TOTEMP.values, LONGLEY_TABLE,
                                 0.015, 0.006, 0.006)
        order = sorted(range(16), key=lambda i: -diags[i].cooks)
        self.assertEqual(order[:2], [4, 15])       # 1951, then 1962
        self.assertEqual(max(range(16), key=lambda i: abs(diags[i].studentized)), 9)  # 1956
        # "approximately the edge of a 15% confidence region" for 1962
        self.assertAlmostEqual(diags[15].level, 0.15, delta=0.02)


class TestReferences(unittest.TestCase):
    def test_statsmodels_and_deletion(self):
        rng = np.random.default_rng(3)
        for n, k in ((13, 4), (30, 2), (9, 3)):
            X0 = rng.normal(size=(n, k))
            y = X0 @ rng.normal(size=k) + rng.normal(size=n)
            y[0] += 4                                       # one outlier
            X, diags = ours(X0, y)
            p = X.shape[1]
            fit = sm.OLS(y, X).fit()
            inf = OLSInfluence(fit)
            np.testing.assert_allclose([d.cooks for d in diags], inf.cooks_distance[0],
                                       rtol=1e-10)
            np.testing.assert_allclose([d.studentized for d in diags],
                                       inf.resid_studentized_internal, rtol=1e-10)
            np.testing.assert_allclose([d.leverage for d in diags], inf.hat_matrix_diag,
                                       rtol=1e-10)
            np.testing.assert_allclose([d.fitted for d in diags], fit.fittedvalues,
                                       rtol=1e-10)
            np.testing.assert_allclose([d.level for d in diags],
                                       f_dist.cdf(inf.cooks_distance[0], p, n - p), atol=1e-9)
            s2 = fit.mse_resid
            for i, d in enumerate(diags):                  # Cook's definition, eq. (4)
                keep = np.arange(n) != i
                delta = sm.OLS(y[keep], X[keep]).fit().params - fit.params
                self.assertAlmostEqual(d.cooks, delta @ X.T @ X @ delta / (p * s2), places=9)
                self.assertEqual(d.influential, d.cooks >= f_dist.ppf(0.5, p, n - p))

    def test_leverage_one_and_perfect_fit(self):
        X = [[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 1.0]]     # last run alone fixes b1
        y = [1.0, 2.0, 3.0, 9.0]
        inv = np.linalg.inv(np.array(X).T @ np.array(X))
        beta = inv @ np.array(X).T @ np.array(y)
        diags = influence.diagnostics(X, y, beta.tolist(), inv.tolist(), 2.0)
        self.assertAlmostEqual(diags[3].leverage, 1.0)
        self.assertIsNone(diags[3].cooks)
        self.assertFalse(diags[3].influential)
        self.assertIsNotNone(diags[0].cooks)
        diags = influence.diagnostics(X, [1.0, 1.0, 1.0, 9.0], [1.0, 8.0], inv.tolist(), 0.0)
        self.assertTrue(all(d.cooks is None for d in diags))


def numeric(k):
    return [doe.Factor(f"X{i + 1}", "", 10.0, 30.0) for i in range(k)]


class TestAnalyses(unittest.TestCase):
    def test_rsm_blocked_with_missing(self):
        d = doe.build_rsm(numeric(3), "ccd", blocks=2, randomize=False)
        coded = [r.coded for r in d.runs]
        blocks = [r.block for r in d.runs]
        rng = np.random.default_rng(5)
        y = [50 + x[0] - x[1] ** 2 + rng.normal() for x in coded]
        y[3] += 25
        y[7] = None
        fit = rsm.fit_model(coded, blocks, y)
        obs = [i for i, v in enumerate(y) if v is not None]
        self.assertEqual(fit.observations, obs)
        X = np.array([[1.0, 1.0 if blocks[i] == 1 else -1.0] + list(quad(coded[i]))[1:]
                      for i in obs])
        inf = OLSInfluence(sm.OLS(np.array([y[i] for i in obs]), X).fit())
        np.testing.assert_allclose([g.cooks for g in fit.diagnostics], inf.cooks_distance[0],
                                   rtol=1e-8)
        self.assertTrue(fit.diagnostics[obs.index(3)].influential)

    def test_screening_and_saturated(self):
        d = doe.build_design(numeric(3), "full", centre_points=3, randomize=False)
        coded = [r.coded for r in d.runs]
        ptype = [r.ptype for r in d.runs]
        y = [10 + 2 * x[0] + x[1] * x[2] + 0.1 * i for i, x in enumerate(coded)]
        a = screening.analyse(coded, ptype, [1] * len(y), y, "main")
        X = np.array([[1.0, 1.0 if p == 0 else 0.0] + list(x) for x, p in zip(coded, ptype)])
        inf = OLSInfluence(sm.OLS(np.array(y), X).fit())
        np.testing.assert_allclose([g.cooks for g in a.diagnostics], inf.cooks_distance[0],
                                   rtol=1e-8, atol=1e-12)
        d = doe.build_design(numeric(2), "full", randomize=False)
        a = screening.analyse([r.coded for r in d.runs], [1] * 4, [1] * 4, [1.0, 2, 4, 3],
                              "2fi")
        self.assertEqual(a.df, 0)
        self.assertIsNone(a.diagnostics)


def quad(x):
    from itertools import combinations
    k = len(x)
    return [1.0] + list(x) + [v * v for v in x] + [x[i] * x[j] for i, j in
                                                  combinations(range(k), 2)]


if __name__ == "__main__":
    unittest.main()


# example_data_rsm.ods, sheet "Design 2": CCCD, 3 factors, 28 runs, Yield (mg/g E)
DESIGN2_Y = [24.28989323, 33.66123277, 39.83253039, 18.65134664, 12.52643918, 24.52400157,
             19.43370304, 18.04370825, 11.58758592, 19.25931338, 16.25878132, 9.714146726,
             13.70971019, 33.93868664, 22.13181921, 51.13647479, 26.76090827, 37.48920235,
             15.37816958, 31.05661974, 11.95990485, 18.31395824, 17.97089328, 11.73768574,
             8.318061785, 12.1663204, 14.31427388, 10.81211525]


class TestOutlierTest(unittest.TestCase):
    def design2(self):
        d = doe.build_rsm([doe.Factor(n, "", 0.0, 1.0) for n in ("Time", "Temp", "EhOH")],
                          "cccd", randomize=False)
        return rsm.fit_model([r.coded for r in d.runs], [1] * 28, DESIGN2_Y)

    def test_design2_validation(self):
        # the user's reference: run 3 -> t* = 4.19, p = 0.0006, p_bonf = 0.017, h = 0.163,
        # D = 0.178 -> "outlier", not influential
        fit = self.design2()
        d = fit.diagnostics[2]
        self.assertAlmostEqual(d.external, 4.19, delta=0.005)
        self.assertAlmostEqual(d.p_outlier, 0.0006, delta=0.00005)
        self.assertAlmostEqual(d.p_bonferroni, 0.017, delta=0.0005)
        self.assertAlmostEqual(d.leverage, 0.163, delta=0.0005)
        self.assertAlmostEqual(d.cooks, 0.178, delta=0.0005)
        self.assertEqual(d.flags(), ["outlier"])
        self.assertFalse(d.influential)
        flags = {i + 1: g.flags() for i, g in enumerate(fit.diagnostics) if g.flags()}
        self.assertEqual(flags, {3: ["outlier"], 9: ["check"]})
        line = influence.verdict("3", d, 28, 10)
        self.assertTrue(line.startswith("Run 3: outlier (t* = 4.19, Bonferroni p = 0.017); "
                                        "not influential (D = 0.18"), line)
        self.assertIn("Action: check the record", line)

    def test_against_statsmodels_and_scipy(self):
        from scipy.stats import t as t_dist
        rng = np.random.default_rng(11)
        for n, k in ((13, 4), (30, 2), (8, 4)):
            X0 = rng.normal(size=(n, k))
            y = X0 @ rng.normal(size=k) + rng.normal(size=n)
            y[1] += 5
            X, diags = ours(X0, y)
            p = X.shape[1]
            inf = OLSInfluence(sm.OLS(y, X).fit())
            np.testing.assert_allclose([d.external for d in diags],
                                       inf.resid_studentized_external, rtol=1e-9)
            pv = 2 * t_dist.sf(np.abs(inf.resid_studentized_external), n - p - 1)
            np.testing.assert_allclose([d.p_outlier for d in diags], pv, rtol=1e-7,
                                       atol=1e-15)
            np.testing.assert_allclose([d.p_bonferroni for d in diags],
                                       np.minimum(1, n * pv), rtol=1e-7, atol=1e-15)
            h = inf.hat_matrix_diag
            self.assertEqual([d.high_leverage for d in diags],
                             list(h > 2 * p / n) if p / n <= 0.5 else [False] * n)

    def test_leverage_rule_and_guards(self):
        self.assertAlmostEqual(influence.leverage_cutoff(28, 10), 20 / 28)
        self.assertIsNone(influence.leverage_cutoff(10, 6))         # p/n > 0.5
        # n - p - 1 = 0: no t*, but r and D are still there
        X = np.column_stack([np.ones(4), [0.0, 1, 2, 3], [0.0, 1, 4, 9]])
        y = np.array([1.0, 2.5, 2.9, 5.0])
        inv = np.linalg.inv(X.T @ X)
        beta = inv @ X.T @ y
        sse = float(((y - X @ beta) ** 2).sum())
        diags = influence.diagnostics(X.tolist(), y.tolist(), beta.tolist(), inv.tolist(), sse)
        self.assertTrue(all(d.external is None and d.cooks is not None for d in diags))
        # r^2 = n - p: one run off a line that the others fit exactly -> t* infinite -> n/a
        X = np.column_stack([np.ones(5), [0.0, 1, 2, 3, 4]])
        y = np.array([0.0, 1, 2, 3, 9])
        inv = np.linalg.inv(X.T @ X)
        beta = inv @ X.T @ y
        sse = float(((y - X @ beta) ** 2).sum())
        diags = influence.diagnostics(X.tolist(), y.tolist(), beta.tolist(), inv.tolist(), sse)
        self.assertAlmostEqual(diags[4].studentized ** 2, 3.0)
        self.assertIsNone(diags[4].external)
        self.assertIsNotNone(diags[0].external)
