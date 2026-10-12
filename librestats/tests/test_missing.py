"""Check librestats.missing against R naniar, scikit-learn and simulation.

Run: ../.venv-dev/bin/python -m unittest discover -s tests -p "test_[mod]*.py"   (from librestats/)
"""

import csv
import os
import random
import statistics
import sys
import unittest

import numpy as np
from sklearn.impute import KNNImputer

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import linalg  # noqa: E402
from librestats import missing as mi  # noqa: E402

HERE = os.path.dirname(__file__)


def airquality():
    with open(os.path.join(HERE, "data", "airquality.csv")) as f:
        reader = csv.reader(f)
        next(reader)
        return [[float(v) if v not in ("", "NA") else None for v in row[1:]] for row in reader]


def with_mcar(rng, n, p=3, rate=0.2, mar=False):
    """Correlated normal data; missingness in columns 1..p-1 either MCAR or
    depending on column 0 (MAR)."""
    cov = np.full((p, p), 0.5) + 0.5 * np.eye(p)
    x = rng.multivariate_normal(np.zeros(p), cov, n)
    rows = x.tolist()
    for i in range(n):
        for j in range(1, p):
            prob = (0.45 if x[i, 0] > 0 else 0.02) if mar else rate
            if rng.random() < prob:
                rows[i][j] = None
    return rows


class TestLinalg(unittest.TestCase):
    def test_against_numpy(self):
        rng = np.random.default_rng(1)
        a = rng.standard_normal((5, 5))
        a = a @ a.T + 5 * np.eye(5)
        b = rng.standard_normal(5)
        L = linalg.cholesky(a.tolist())
        np.testing.assert_allclose(L, np.linalg.cholesky(a), rtol=1e-12)
        np.testing.assert_allclose(linalg.cho_solve(L, b.tolist()), np.linalg.solve(a, b),
                                   rtol=1e-10)
        np.testing.assert_allclose(linalg.inverse_spd(a.tolist()), np.linalg.inv(a), rtol=1e-10)
        self.assertAlmostEqual(linalg.quad_form_inv(a.tolist(), b.tolist()),
                               b @ np.linalg.solve(a, b), places=10)

    def test_singular(self):
        with self.assertRaises(linalg.SingularMatrixError):
            linalg.cholesky([[1.0, 2.0], [2.0, 4.0]])


class TestPatternsAndEM(unittest.TestCase):
    def test_patterns(self):
        rows = [[1, None], [2, 3], [None, 4], [5, 6], [7, None]]
        pts = mi.patterns(rows)
        self.assertEqual([(pt.observed, pt.count) for pt in pts],
                         [((True, True), 2), ((True, False), 2), ((False, True), 1)])

    def test_em_complete_data_equals_ml_estimates(self):
        x = np.random.default_rng(2).standard_normal((50, 3))
        em = mi.em_mvnorm(x.tolist())
        np.testing.assert_allclose(em.mean, x.mean(0), rtol=1e-12)
        np.testing.assert_allclose(em.cov, np.cov(x, rowvar=False, ddof=0), rtol=1e-10)


class TestLittle(unittest.TestCase):
    def test_airquality_matches_naniar(self):
        # naniar tests/testthat/test-mcar-test.R: mcar_test(airquality).
        # naniar's EM (norm::em.norm) stops early; fully converged ML gives
        # 35.10618 (confirmed by direct optimisation, tools/check_little_ml.py),
        # so allow 1e-5 relative difference.
        res = mi.little_mcar(airquality())
        self.assertAlmostEqual(res.statistic, 35.1061288689702, delta=35.1061 * 1e-5)
        self.assertEqual(res.df, 14)
        self.assertAlmostEqual(res.p_value, 0.00141778113856683, delta=1e-7)
        self.assertEqual(res.n_patterns, 4)

    def test_all_missing_row_is_excluded(self):
        rows = airquality() + [[None] * 6]
        res = mi.little_mcar(rows)
        self.assertEqual(res.excluded_rows, 1)
        self.assertAlmostEqual(res.statistic, 35.10613, places=4)

    def test_null_and_power(self):
        rng = np.random.default_rng(3)
        reps = 200
        null = sum(mi.little_mcar(with_mcar(rng, 150)).p_value < 0.05 for _ in range(reps))
        self.assertAlmostEqual(null / reps, 0.05, delta=0.035)
        power = sum(mi.little_mcar(with_mcar(rng, 150, mar=True)).p_value < 0.05
                    for _ in range(50))
        self.assertGreater(power / 50, 0.9)

    def test_no_missing_raises(self):
        with self.assertRaises(mi.MissingDataError):
            mi.little_mcar([[1.0, 2.0], [2.0, 1.0], [3.0, 5.0]])


class TestSimpleAndKNN(unittest.TestCase):
    def test_mean_median(self):
        rows = [[1.0, 10.0], [None, 20.0], [3.0, None], [8.0, 30.0]]
        mean = mi.impute_simple(rows, "mean")
        self.assertEqual(mean.data[1][0], 4.0)
        self.assertEqual(mean.data[2][1], 20.0)
        self.assertEqual(sorted(mean.imputed), [(1, 0), (2, 1)])
        self.assertEqual(mi.impute_simple(rows, "median").data[1][0], 3.0)

    def _sklearn(self, x, k):
        return KNNImputer(n_neighbors=k).fit_transform(x)

    def test_knn_matches_sklearn_raw(self):
        rows = with_mcar(np.random.default_rng(4), 80, p=4)
        x = np.array([[np.nan if v is None else v for v in r] for r in rows])
        for k in (1, 3, 5):
            ours = mi.impute_knn(rows, k=k, standardize=False).data
            np.testing.assert_allclose(ours, self._sklearn(x, k), rtol=1e-12)

    def test_knn_matches_sklearn_standardized(self):
        rows = with_mcar(np.random.default_rng(5), 80, p=4)
        x = np.array([[np.nan if v is None else v for v in r] for r in rows])
        x[:, 2] *= 1000  # very different scales
        rows = [[None if np.isnan(v) else float(v) for v in r] for r in x]
        mu, sd = np.nanmean(x, 0), np.nanstd(x, 0, ddof=1)
        ref = self._sklearn((x - mu) / sd, 5) * sd + mu
        np.testing.assert_allclose(mi.impute_knn(rows, k=5).data, ref, rtol=1e-10)


class TestMice(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # y depends on x (rho = 0.6); y is missing more often when x is high (MAR)
        rng = np.random.default_rng(6)
        n = 600
        x = rng.standard_normal(n)
        y = 0.6 * x + 0.8 * rng.standard_normal(n)
        cls.true_mean = 0.0
        rows = []
        for xi, yi in zip(x, y):
            rows.append([float(xi), None if rng.random() < (0.6 if xi > 0 else 0.1) else float(yi)])
        cls.rows = rows
        cls.cc_mean = statistics.fmean(r[1] for r in rows if r[1] is not None)

    def check(self, method):
        res = mi.mice(self.rows, m=5, iterations=10, method=method, seed=42)
        self.assertEqual(len(res.datasets), 5)
        means = [statistics.fmean(d[1] for d in ds) for ds in res.datasets]
        sds = [statistics.stdev(d[1] for d in ds) for ds in res.datasets]
        cors = [np.corrcoef(np.array(ds).T)[0, 1] for ds in res.datasets]
        # complete-case mean is biased downwards; MICE should remove most of that
        self.assertLess(self.cc_mean, -0.15)
        self.assertLess(abs(statistics.fmean(means)), 0.08)
        self.assertAlmostEqual(statistics.fmean(sds), 1.0, delta=0.1)
        self.assertAlmostEqual(statistics.fmean(cors), 0.6, delta=0.08)
        self.assertGreater(statistics.stdev(means), 0)  # imputations differ
        return res

    def test_pmm(self):
        res = self.check("pmm")
        observed = {r[1] for r in self.rows if r[1] is not None}
        self.assertTrue(all(res.datasets[0][i][j] in observed for i, j in res.imputed))

    def test_norm(self):
        self.check("norm")

    def test_seed_reproducible(self):
        a = mi.mice(self.rows[:100], m=2, iterations=3, seed=7).datasets
        b = mi.mice(self.rows[:100], m=2, iterations=3, seed=7).datasets
        self.assertEqual(a, b)

    def test_observed_values_unchanged(self):
        res = mi.mice(self.rows[:100], m=2, iterations=3, seed=1)
        for ds in res.datasets:
            for r_orig, r_new in zip(self.rows[:100], ds):
                for v_orig, v_new in zip(r_orig, r_new):
                    if v_orig is not None:
                        self.assertEqual(v_orig, v_new)


if __name__ == "__main__":
    unittest.main()
