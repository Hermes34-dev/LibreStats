"""Check librestats.histogram by hand-counted examples and against numpy.histogram."""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import histogram as h  # noqa: E402

RNG = np.random.default_rng(2026)
ZERO_TO_TEN = [float(i) for i in range(11)]


class TestHandCounts(unittest.TestCase):
    def test_number_of_bins_with_values_on_edges(self):
        # edges 0, 2, 4, 6, 8, 10: 2, 4, 6, 8 start the next bin, 10 stays in the last
        res = h.histogram(ZERO_TO_TEN, "count", k=5)
        self.assertEqual(res.edges, [0, 2, 4, 6, 8, 10])
        self.assertEqual(res.counts, [2, 2, 2, 2, 3])
        self.assertEqual(res.labels(), ["0–2", "2–4", "4–6", "6–8", "8–10"])
        self.assertEqual(res.width, 2)

    def test_bin_width(self):
        # k = ceil(10 / 3) = 4, edges 0, 3, 6, 9, 12
        res = h.histogram(ZERO_TO_TEN, "width", width=3)
        self.assertEqual(res.edges, [0, 3, 6, 9, 12])
        self.assertEqual(res.counts, [3, 3, 3, 2])

    def test_decimal_width_edges(self):
        # 0.3 / 0.1 is 2.9999999999999996 in floating point; 0.3 must still start bin 4
        xs = [0.0, 0.1, 0.2, 0.3, 0.3, 0.7, 1.0]
        res = h.histogram(xs, "width", width=0.1)
        self.assertEqual(len(res.counts), 10)
        self.assertEqual(res.counts, [1, 1, 1, 2, 0, 0, 0, 1, 0, 1])
        self.assertEqual(res.labels()[3], "0.3–0.4")

    def test_sturges(self):
        self.assertEqual([h.sturges(n) for n in (2, 8, 11, 16, 100)], [2, 4, 5, 5, 8])
        res = h.histogram(ZERO_TO_TEN)
        self.assertEqual(len(res.counts), 5)   # ceil(log2 11) + 1
        self.assertEqual(res.counts, [2, 2, 2, 2, 3])  # edges 0, 2, ..., 10 as with k = 5
        self.assertIn("Sturges", res.rule)

    def test_user_range(self):
        # only 2..8 counted (7 values), edges 2, 4, 6, 8
        res = h.histogram(ZERO_TO_TEN, "count", k=3, lo=2, hi=8)
        self.assertEqual(res.counts, [2, 2, 3])
        self.assertEqual(res.excluded, 4)
        self.assertEqual(res.n, 7)

    def test_range_wider_than_data(self):
        res = h.histogram([1.0, 2.0, 3.0], "width", width=5, lo=0, hi=10)
        self.assertEqual(res.edges, [0, 5, 10])
        self.assertEqual(res.counts, [3, 0])

    def test_only_min_given(self):
        res = h.histogram(ZERO_TO_TEN, "count", k=2, lo=5)
        self.assertEqual(res.edges, [5, 7.5, 10])
        self.assertEqual(res.counts, [3, 3])
        self.assertEqual(res.excluded, 5)

    def test_percent(self):
        res = h.histogram(ZERO_TO_TEN, "count", k=5)
        self.assertAlmostEqual(sum(res.percents()), 100.0)
        self.assertAlmostEqual(res.percents()[4], 300 / 11)

    def test_negative_values(self):
        res = h.histogram([-3.0, -1.0, 0.0, 1.0], "count", k=2)
        self.assertEqual(res.edges, [-3, -1, 1])
        self.assertEqual(res.counts, [1, 3])
        self.assertEqual(res.labels(), ["-3–-1", "-1–1"])


class TestEdgeCases(unittest.TestCase):
    def test_empty(self):
        with self.assertRaises(h.HistogramError):
            h.histogram([])

    def test_single_value(self):
        res = h.histogram([5.0])
        self.assertEqual((res.edges, res.counts), ([5.0, 5.0], [1]))
        self.assertEqual(res.labels(), ["5–5"])
        self.assertTrue(res.notes)

    def test_identical_values(self):
        for bins, k in (("auto", None), ("count", 7)):
            res = h.histogram([2.5] * 4, bins, k=k)
            self.assertEqual(res.counts, [4])
        res = h.histogram([2.5] * 4, "width", width=0.5)
        self.assertEqual((res.edges, res.counts), ([2.5, 3.0], [4]))

    def test_invalid_settings(self):
        for kwargs in ({"bins": "count", "k": 0}, {"bins": "count", "k": 2.5},
                       {"bins": "width", "width": 0}, {"bins": "width", "width": -1},
                       {"bins": "width", "width": 1e-6},   # too many bins
                       {"lo": 5, "hi": 5}, {"lo": 6, "hi": 2},
                       {"lo": 20}):                         # nothing in range
            with self.assertRaises(h.HistogramError, msg=str(kwargs)):
                h.histogram(ZERO_TO_TEN, **kwargs)

    def test_numeric_values_flattens_and_skips_text(self):
        data = ((1.0, "a", ""), (2.0, "  ", 3.0), ("NA", 4.0, float("nan")))
        values, text = h.numeric_values(data)
        self.assertEqual(values, [1.0, 2.0, 3.0, 4.0])
        self.assertEqual(text, 2)


class TestAgainstNumpy(unittest.TestCase):
    """numpy.histogram uses the same convention: [a, b), last bin closed."""

    def test_number_of_bins(self):
        for n, k in ((10, 3), (57, 7), (500, 12), (1000, 1)):
            x = RNG.normal(10, 3, n)
            res = h.histogram(list(x), "count", k=k)
            counts, edges = np.histogram(x, bins=k)
            self.assertEqual(res.counts, list(counts))
            np.testing.assert_allclose(res.edges, edges, rtol=1e-10)

    def test_integer_data_on_edges(self):
        x = RNG.integers(0, 20, 300).astype(float)
        res = h.histogram(list(x), "count", k=10, lo=0, hi=20)
        counts, _ = np.histogram(x, bins=10, range=(0, 20))
        self.assertEqual(res.counts, list(counts))

    def test_bin_width(self):
        x = RNG.uniform(0, 7.3, 200)
        res = h.histogram(list(x), "width", width=0.5, lo=0)
        counts, _ = np.histogram(x, bins=np.arange(0, 7.5 + 1e-9, 0.5))
        self.assertEqual(res.counts, list(counts))

    def test_sturges_matches_numpy(self):
        for n in (5, 30, 129, 1000):
            x = RNG.normal(0, 1, n)
            self.assertEqual(len(h.histogram(list(x)).counts),
                             len(np.histogram_bin_edges(x, bins="sturges")) - 1)


if __name__ == "__main__":
    unittest.main()
