"""Check librestats.coding against pandas, scikit-learn and brute force."""

import os
import sys
import unittest

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import coding as cd  # noqa: E402

RNG = np.random.default_rng(99)


def brute_kmeans_ss(x, k):
    """Plain O(k n^2) dynamic programme for the optimal 1-D k-means SS."""
    x = sorted(x)
    n = len(x)

    def ss(j, i):
        seg = x[j:i + 1]
        m = sum(seg) / len(seg)
        return sum((v - m) ** 2 for v in seg)

    D = [[np.inf] * n for _ in range(k)]
    for i in range(n):
        D[0][i] = ss(0, i)
    for m in range(1, k):
        for i in range(m, n):
            D[m][i] = min(D[m - 1][j - 1] + ss(j, i) for j in range(m, i + 1))
    return D[k - 1][n - 1]


class TestDummiesAndPresence(unittest.TestCase):
    def test_dummies_match_pandas(self):
        vals = ["b", "a", None, "c", "a", "", "b"]
        res = cd.dummies(vals)
        self.assertEqual(res.names, ["a", "b", "c"])
        ref = pd.get_dummies(pd.Series(vals).replace("", None), dtype=float)
        for name, col in zip(res.names, res.columns):
            for v, got, want in zip(vals, col, ref[name]):
                self.assertEqual(got, None if cd.is_missing(v) else want)

    def test_drop_first_and_numbers(self):
        res = cd.dummies([2.0, 1.0, 3.0, 1.0], drop_first=True)
        self.assertEqual(res.names, ["2", "3"])
        self.assertEqual(res.columns[0], [1.0, 0.0, 0.0, 0.0])
        self.assertEqual(dict(res.params)["Reference (no column)"], "1")

    def test_presence(self):
        res = cd.presence_absence([0.0, 3.0, None, "x", 0.5], threshold=0.5)
        self.assertEqual(res.columns[0], [0.0, 1.0, None, None, 0.0])
        res = cd.presence_absence([0.0, None], empty_is_absent=True)
        self.assertEqual(res.columns[0], [0.0, 0.0])


class TestBins(unittest.TestCase):
    def setUp(self):
        self.x = RNG.gamma(2.0, 3.0, 200).round(2).tolist()

    def test_equal_width_matches_pandas_cut(self):
        for k in (2, 4, 7):
            codes = cd.bins_equal_width(self.x + [None], k).columns[0]
            ref = pd.cut(self.x, k, labels=False) + 1
            self.assertEqual(codes[:-1], ref.tolist())
            self.assertIsNone(codes[-1])

    def test_quantile_matches_pandas_qcut(self):
        for k in (2, 4, 5, 10):
            codes = cd.bins_quantile(self.x, k).columns[0]
            ref = pd.qcut(self.x, k, labels=False, duplicates="drop") + 1
            self.assertEqual(codes, ref.tolist())

    def test_quantile_ties(self):
        res = cd.bins_quantile([1.0] * 10 + [2.0, 3.0], 4)
        self.assertTrue(res.notes)

    def test_kmeans_optimal(self):
        x = RNG.normal(0, 1, 40).tolist()
        for k in (2, 3, 5):
            _, _, ss = cd.kmeans_1d(x, k)
            self.assertAlmostEqual(ss, brute_kmeans_ss(x, k), places=9)

    def test_kmeans_not_worse_than_sklearn(self):
        x = np.concatenate([RNG.normal(m, 0.7, 80) for m in (0, 4, 9, 15)])
        for k in (3, 4, 6):
            _, _, ss = cd.kmeans_1d(x.tolist(), k)
            km = KMeans(k, n_init=20, random_state=0).fit(x.reshape(-1, 1))
            self.assertLessEqual(ss, km.inertia_ * (1 + 1e-9))

    def test_kmeans_bins_are_contiguous(self):
        res = cd.bins_kmeans(self.x, 4)
        codes = res.columns[0]
        pairs = sorted(zip(self.x, codes))
        self.assertEqual([c for _, c in pairs], sorted(c for _, c in pairs))
        self.assertEqual(sum(row[3] for row in res.table[1:]), len(self.x))


class TestRecode(unittest.TestCase):
    def test_rules(self):
        rules = cd.parse_rules("""
            # merge Likert categories
            1; 2 -> low
            3 -> mid
            4;5 -> high
        """)
        res = cd.recode([1.0, 3.0, 5.0, None, 2.0, 6.0], rules)
        self.assertEqual(res.columns[0], ["low", "mid", "high", None, "low", 6.0])
        self.assertIn("6", res.notes[0])

    def test_text_and_default(self):
        rules = cd.parse_rules("yes; Y -> 1\nno -> 0\n* -> other")
        res = cd.recode(["yes", "no", "maybe", "Y", ""], rules)
        self.assertEqual(res.columns[0], [1.0, 0.0, "other", 1.0, None])

    def test_template_and_errors(self):
        self.assertEqual(cd.rules_template(["b", 2.0, "a", 2.0, 10.0]),
                         "2 -> 2\n10 -> 10\na -> a\nb -> b")
        with self.assertRaises(cd.CodingError):
            cd.parse_rules("1 = low")


if __name__ == "__main__":
    unittest.main()
