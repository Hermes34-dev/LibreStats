"""Check librestats.compositional against scikit-bio and the defining properties."""

import os
import sys
import unittest

import numpy as np
from skbio.stats import composition as skc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import compositional as c  # noqa: E402

RNG = np.random.default_rng(11)
X = RNG.dirichlet([5, 2, 8, 1, 3], 30) * 100  # percentages summing to 100
NAMES = ["C16:0", "C18:0", "C18:1", "C18:2", "C18:3"]


class TestLogRatios(unittest.TestCase):
    def test_clr_matches_skbio(self):
        res = c.log_ratios(X.tolist(), NAMES, "clr")
        np.testing.assert_allclose(res.rows, skc.clr(X), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(np.sum(res.rows, 1), 0, atol=1e-12)

    def test_ilr_matches_skbio_with_pivot_basis(self):
        res = c.log_ratios(X.tolist(), NAMES, "ilr")
        basis = np.array(c.pivot_basis(5))
        np.testing.assert_allclose(basis @ basis.T, np.eye(4), atol=1e-12)  # orthonormal
        np.testing.assert_allclose(basis.sum(1), 0, atol=1e-12)             # clr contrasts
        np.testing.assert_allclose(res.rows, skc.ilr(X, basis=basis), rtol=1e-10, atol=1e-12)

    def test_ilr_is_isometric(self):
        z = np.array(c.log_ratios(X.tolist(), NAMES, "ilr").rows)
        k = np.array(c.log_ratios(X.tolist(), NAMES, "clr").rows)
        for i, j in ((0, 1), (3, 17), (5, 29)):
            self.assertAlmostEqual(np.linalg.norm(z[i] - z[j]), np.linalg.norm(k[i] - k[j]),
                                   places=10)

    def test_scale_invariance(self):
        a = c.log_ratios(X.tolist(), NAMES, "ilr").rows
        b = c.log_ratios((X / 100).tolist(), NAMES, "ilr").rows
        np.testing.assert_allclose(a, b, atol=1e-12)

    def test_missing_rows_skipped(self):
        rows = X[:3].tolist()
        rows[1][2] = None
        res = c.log_ratios(rows, NAMES, "clr")
        self.assertIsNone(res.rows[1])
        self.assertIsNotNone(res.rows[0])

    def test_zeros(self):
        rows = [[0.0, 30.0, 70.0], [10.0, 20.0, 70.0], [5.0, 0.0, 95.0]]
        with self.assertRaises(c.CompositionError):
            c.log_ratios(rows, ["a", "b", "c"], "clr")
        new, deltas, replaced = c.replace_zeros(rows)
        self.assertEqual(replaced, 2)
        self.assertEqual(deltas[:2], [0.65 * 5.0, 0.65 * 20.0])
        self.assertAlmostEqual(sum(new[0]), 100.0)               # total kept
        self.assertAlmostEqual(new[0][1] / new[0][2], 30 / 70)   # ratios kept
        self.assertEqual(new[1], rows[1])                        # rows without zeros unchanged
        res = c.log_ratios(rows, ["a", "b", "c"], "clr", zeros="replace")
        self.assertEqual(dict(res.params)["Zeros replaced"], 2)

    def test_negative_raises(self):
        with self.assertRaises(c.CompositionError):
            c.log_ratios([[1.0, -2.0], [1.0, 2.0]], ["a", "b"])


if __name__ == "__main__":
    unittest.main()
