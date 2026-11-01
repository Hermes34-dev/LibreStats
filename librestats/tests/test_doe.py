"""Check librestats.doe: fractions against NIST and Montgomery, Plackett-Burman
orthogonality, blocking against Montgomery's Table 7.8, and the run layout."""

import math
import os
import sys
import unittest
from collections import Counter
from itertools import combinations

import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src", "pythonpath"))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from librestats import doe  # noqa: E402
from librestats.ff_table import DESIGNS  # noqa: E402
import make_ff_table  # noqa: E402


def wlp(fraction):
    counts = Counter(bin(w).count("1") for w in fraction.words)
    return tuple(counts.get(n, 0) for n in range(3, fraction.k + 1))


def wlp_of(k, generators):
    """Word-length pattern of a design given as {added letter: product}."""
    words = doe.subgroup([doe._mask(added + product) for added, product in generators.items()])
    counts = Counter(bin(w).count("1") for w in words)
    return tuple(counts.get(n, 0) for n in range(3, k + 1))


def factors(k, low=-1.0, high=1.0):
    return [doe.Factor(f"X{i + 1}", "", low, high) for i in range(k)]


class TestFractions(unittest.TestCase):
    def test_resolutions_match_nist(self):
        self.assertEqual({key: res for key, (res, _, _) in DESIGNS.items()}, make_ff_table.NIST)

    def test_table_is_minimum_aberration(self):
        # a fresh exhaustive search gives the same word-length pattern
        for (k, runs), (_, _, pattern) in DESIGNS.items():
            if runs <= 64 and k <= 9:
                found, _ = make_ff_table.search(k, runs.bit_length() - 1)
                self.assertEqual(tuple(found[3:]), pattern, (k, runs))

    def test_designs_are_regular_fractions(self):
        for (k, runs), (res, _, pattern) in DESIGNS.items():
            f = doe.fractional_factorial(k, runs)
            x = np.array(f.matrix())
            self.assertEqual(x.shape, (runs, k))
            self.assertEqual(len({tuple(r) for r in x}), runs)            # no repeated runs
            self.assertTrue(np.all(x.sum(axis=0) == 0))                    # balanced
            self.assertTrue(np.array_equal(x.T @ x, runs * np.eye(k)))    # orthogonal
            for w in f.words:                                              # I = word on every run
                cols = [i for i in range(k) if w >> i & 1]
                self.assertTrue(np.all(np.prod(x[:, cols], axis=1) == 1))
            self.assertEqual(len(f.words), 2 ** (k - (runs.bit_length() - 1)) - 1)
            self.assertEqual(min(bin(w).count("1") for w in f.words), res)
            self.assertEqual(wlp(f), pattern)

    def test_nist_2_7_3_defining_relation(self):
        # NIST e-Handbook 2to7m3: I = 1235 = 2346 = 1347 = 1456 = 2457 = 1267 = 3567
        nist = ["1235", "2346", "1347", "1456", "2457", "1267", "3567"]
        counts = Counter(len(w) for w in nist)
        self.assertEqual(wlp(doe.fractional_factorial(7, 16)), (0, counts[4], 0, 0, 0))

    def test_montgomery_generators_same_aberration(self):
        # Montgomery, Design and Analysis of Experiments, Table 8.14
        published = {
            (5, 8): {"D": "AB", "E": "AC"},
            (6, 16): {"E": "ABC", "F": "BCD"},
            (7, 32): {"F": "ABCD", "G": "ABDE"},
            (7, 16): {"E": "ABC", "F": "BCD", "G": "ACD"},
            (8, 16): {"E": "BCD", "F": "ACD", "G": "ABC", "H": "ABD"},
            (9, 32): {"F": "BCDE", "G": "ACDE", "H": "ABDE", "J": "ABCE"},
            (10, 32): {"F": "ABCD", "G": "ABCE", "H": "ABDE", "J": "ACDE", "K": "BCDE"},
        }
        for (k, runs), gens in published.items():
            self.assertEqual(wlp(doe.fractional_factorial(k, runs)), wlp_of(k, gens), (k, runs))

    def test_alias_chains_2_4_1(self):
        d = doe.build_design(factors(4), "fractional", 8, randomize=False)
        self.assertEqual(d.generator_text(), ["D = ABC"])
        self.assertEqual(d.defining_relation(), "I = ABCD")
        self.assertEqual(d.alias_chains(), ["A + BCD", "B + ACD", "C + ABD", "D + ABC",
                                            "AB + CD", "AC + BD", "AD + BC"])

    def test_default_runs(self):
        self.assertEqual(doe.fractional_factorial(5).runs, 16)   # V rather than III
        self.assertEqual(doe.fractional_factorial(7).runs, 16)   # smallest IV
        self.assertEqual(doe.fractional_factorial(3).runs, 4)    # only III exists
        with self.assertRaisesRegex(doe.DesignError, "8 runs .resolution III., 16 runs"):
            doe.fractional_factorial(6, 64)
        with self.assertRaises(doe.DesignError):
            doe.fractional_factorial(12)


class TestPlackettBurman(unittest.TestCase):
    def test_orthogonal_for_all_sizes(self):
        for n in doe.PB_ROWS:
            x = np.array(doe.plackett_burman(n - 1, n))
            self.assertEqual(x.shape, (n, n - 1))
            self.assertTrue(np.array_equal(x.T @ x, n * np.eye(n - 1)), n)

    def test_size_choice(self):
        self.assertEqual(len(doe.plackett_burman(7)), 8)
        self.assertEqual(len(doe.plackett_burman(8)), 12)
        self.assertEqual(len(doe.plackett_burman(11, 20)), 20)
        self.assertEqual(doe.plackett_burman(11)[0], [1, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1])
        with self.assertRaises(doe.DesignError):
            doe.plackett_burman(24)
        with self.assertRaises(doe.DesignError):
            doe.plackett_burman(11, 8)


class TestBlocks(unittest.TestCase):
    def pattern(self, k, blocks):
        d = doe.build_design(factors(k), "full", blocks=blocks, randomize=False)
        return sorted(bin(e).count("1") for e in doe.subgroup(d.block_effects))

    def test_full_factorial_as_montgomery_table_7_8(self):
        # orders of the effects confounded with blocks by Montgomery's generators
        table = {(3, 2): ["ABC"], (3, 4): ["AB", "AC"], (4, 2): ["ABCD"],
                 (4, 4): ["ABC", "ACD"], (4, 8): ["AB", "BC", "CD"], (5, 2): ["ABCDE"],
                 (5, 4): ["ABC", "CDE"], (5, 8): ["ABE", "BCE", "CDE"],
                 (6, 2): ["ABCDEF"], (6, 4): ["ABCF", "CDEF"], (6, 8): ["ABEF", "ABCD", "ACE"]}
        for (k, blocks), gens in table.items():
            ref = sorted(bin(e).count("1") for e in doe.subgroup([doe._mask(g) for g in gens]))
            self.assertEqual(self.pattern(k, blocks), ref, (k, blocks))

    def test_blocks_are_orthogonal_to_main_effects(self):
        for kind, k, runs, blocks in (("full", 5, None, 8), ("fractional", 7, 16, 2),
                                      ("fractional", 6, 16, 4), ("fractional", 11, 128, 8)):
            d = doe.build_design(factors(k), kind, runs, blocks=blocks, randomize=False)
            self.assertEqual(sorted(Counter(r.block for r in d.runs).values()),
                             [d.base_runs // blocks] * blocks)
            for b in range(1, blocks + 1):
                x = np.array([r.coded for r in d.runs if r.block == b])
                self.assertTrue(np.all(x.sum(axis=0) == 0), (kind, k, b))
            for e in doe.subgroup(d.block_effects):
                self.assertTrue(all(bin(a).count("1") > 1 for a in d.fraction.aliases(e)))

    def test_fraction_blocks_avoid_two_factor_interactions_when_possible(self):
        d = doe.build_design(factors(6), "fractional", 32, blocks=2, randomize=False)
        self.assertTrue(all(bin(a).count("1") >= 3 for e in d.block_effects
                            for a in d.fraction.aliases(e)))


class TestBuild(unittest.TestCase):
    def test_runs_centre_points_replicates(self):
        d = doe.build_design(factors(3, 20.0, 80.0), "full", centre_points=3, replicates=2,
                             blocks=2, seed=7)
        self.assertEqual(len(d.runs), 2 * 8 + 2 * 3)
        self.assertEqual([r.order for r in d.runs], list(range(1, 23)))
        self.assertEqual(sorted(r.std for r in d.runs), list(range(1, 23)))
        self.assertEqual([r.block for r in d.runs], [1] * 11 + [2] * 11)   # replicate = block
        self.assertEqual(sum(r.centre for r in d.runs if r.block == 2), 3)
        self.assertEqual({tuple(r.coded) for r in d.runs if r.centre}, {(0, 0, 0)})
        f = d.factors[0]
        self.assertEqual((f.actual(-1), f.actual(0), f.actual(1)), (20.0, 50.0, 80.0))

    def test_randomisation_is_reproducible_and_within_blocks(self):
        a = doe.build_design(factors(4), "full", blocks=4, seed=42)
        b = doe.build_design(factors(4), "full", blocks=4, seed=42)
        c = doe.build_design(factors(4), "full", blocks=4, seed=43)
        self.assertEqual([r.std for r in a.runs], [r.std for r in b.runs])
        self.assertNotEqual([r.std for r in a.runs], [r.std for r in c.runs])
        self.assertEqual([r.block for r in a.runs], sorted(r.block for r in a.runs))
        d = doe.build_design(factors(4), "full", randomize=False)
        self.assertIsNone(d.seed)
        self.assertEqual([r.std for r in d.runs], list(range(1, 17)))
        self.assertEqual(d.runs[1].coded, [1, -1, -1, -1])            # Yates order
        self.assertIsNotNone(doe.build_design(factors(3), "full").seed)  # random seed reported

    def test_text_levels(self):
        fs = factors(2) + [doe.Factor("Catalyst", "", "Pd", "Pt")]
        d = doe.build_design(fs, "full", randomize=False)
        self.assertEqual([fs[2].actual(r.coded[2]) for r in d.runs], ["Pd"] * 4 + ["Pt"] * 4)
        with self.assertRaisesRegex(doe.DesignError, "Catalyst has text levels"):
            doe.build_design(fs, "full", centre_points=1)

    def test_block_errors(self):
        with self.assertRaisesRegex(doe.DesignError, "multiple of the replicates"):
            doe.build_design(factors(3), "full", replicates=2, blocks=3)
        with self.assertRaisesRegex(doe.DesignError, "1, 2, 4 or 8"):
            doe.build_design(factors(5), "full", blocks=16)
        with self.assertRaisesRegex(doe.DesignError, "blocked by replicate"):
            doe.build_design(factors(5), "pb", blocks=2)
        with self.assertRaisesRegex(doe.DesignError, "main effect"):
            doe.build_design(factors(7), "fractional", 8, blocks=2)    # saturated
        d = doe.build_design(factors(5), "pb", replicates=2, blocks=2, randomize=False)
        self.assertEqual([r.block for r in d.runs], [1] * 8 + [2] * 8)

    def test_pb_notes(self):
        self.assertIn("partial aliasing", doe.build_design(factors(5), "pb", 12).notes[0])
        self.assertIn("regular", doe.build_design(factors(5), "pb", 8).notes[0])


def quad(x):
    k = len(x)
    return [1.0] + list(x) + [v * v for v in x] + [x[i] * x[j] for i, j in combinations(range(k), 2)]


def prediction_variance(points, x):
    X = np.array([quad(p) for p in points], float)
    v = np.array(quad(x))
    return v @ np.linalg.solve(X.T @ X, v)


def numeric(k):
    return [doe.Factor(f"X{i + 1}", "", 10.0, 30.0) for i in range(k)]


class TestCentralComposite(unittest.TestCase):
    def test_uniform_precision_table(self):
        # Box & Hunter (1957), as tabulated in Montgomery (Table 11.x): (k, cube runs, centre)
        for k, cube, centre in ((2, 4, 5), (3, 8, 6), (4, 16, 7), (5, 16, 6), (5, 32, 10),
                                (6, 32, 9), (6, 64, 15), (7, 64, 14), (7, 128, 21)):
            self.assertEqual(doe.uniform_precision_centre(k, cube), centre, (k, cube))
            d = doe.build_rsm(numeric(k), "ccd", fraction=cube < 2 ** k, randomize=False)
            self.assertEqual(len(d.runs), cube + 2 * k + centre)
            self.assertAlmostEqual(d.alpha, cube ** 0.25)

    def test_rotatable_by_prediction_variance(self):
        for k in (2, 3, 5):
            d = doe.build_rsm(numeric(k), "ccd", randomize=False)
            pts = [r.coded for r in d.runs]
            rng = np.random.default_rng(k)
            for radius in (0.5, 1.3):
                dirs = rng.normal(size=(5, k))
                v = [prediction_variance(pts, radius * u / np.linalg.norm(u)) for u in dirs]
                self.assertLess(max(v) - min(v), 1e-9 * max(v))
            self.assertEqual(dict(d.properties)["Rotatable"], "yes")
        d = doe.build_rsm(numeric(3), "ccd", alpha="face", randomize=False)
        pts = [r.coded for r in d.runs]
        self.assertNotAlmostEqual(prediction_variance(pts, [1, 0, 0]),
                                  prediction_variance(pts, [3 ** -0.5] * 3))
        self.assertEqual(dict(d.properties)["Rotatable"], "no")
        self.assertEqual({abs(v) for r in d.runs for v in r.coded}, {0, 1})

    def test_axial_points_and_layout(self):
        d = doe.build_rsm(numeric(3), "ccd", randomize=False)
        axial = [r.coded for r in d.runs if r.ptype == doe.AXIAL]
        self.assertEqual(len(axial), 6)
        self.assertTrue(all(sorted(map(abs, a)) == [0, 0, d.alpha] for a in axial))
        self.assertEqual([r.ptype for r in d.runs], [1] * 8 + [-1] * 6 + [0] * 6)
        f = d.factors[0]
        self.assertAlmostEqual(f.actual(-d.alpha), 20 - 10 * d.alpha)

    def test_inscribed(self):
        d = doe.build_rsm(numeric(2), "ccd", inscribed=True, randomize=False)
        f = d.factors[0]
        self.assertEqual(d.scale, d.alpha)
        self.assertAlmostEqual(f.actual(-d.alpha, d.scale), 10.0)
        self.assertAlmostEqual(f.actual(d.alpha, d.scale), 30.0)
        self.assertAlmostEqual(f.actual(1, d.scale), 20 + 10 / 2 ** 0.5)
        self.assertFalse(any("outside" in n for n in d.notes))

    def test_orthogonal_blocking_textbook_designs(self):
        # Montgomery / Box & Hunter: k=2 in 2 blocks with 3 + 3 centre points: alpha = sqrt 2,
        # rotatable and orthogonal; k=3 in 3 blocks, 2 centre points each: alpha = 1.633;
        # k=4 in 3 blocks, 2 each: alpha = 2, rotatable and orthogonal
        for k, blocks, cf, ca, alpha, rotatable in ((2, 2, 3, 3, 2 ** 0.5, "yes"),
                                                    (3, 3, 2, 2, 1.6330, "no"),
                                                    (4, 3, 2, 2, 2.0, "yes")):
            d = doe.build_rsm(numeric(k), "ccd", alpha="orthogonal", centre_points=cf,
                              axial_centre_points=ca, blocks=blocks, randomize=False)
            self.assertAlmostEqual(d.alpha, alpha, places=4)
            props = dict(d.properties)
            self.assertEqual(props["Blocks orthogonal to the quadratic model"], "yes")
            self.assertEqual(props["Rotatable"], rotatable)
            self.assertTrue(block_dummies_orthogonal(d))
            self.assertEqual(d.block_text(), [] if blocks == 2 else
                             ["ABC" if k == 3 else "ABCD"])

    def test_rotatable_alpha_in_blocks_may_not_be_orthogonal(self):
        d = doe.build_rsm(numeric(3), "ccd", blocks=2, randomize=False)
        self.assertEqual((d.centre_points, d.axial_centre_points), (2, 2))
        self.assertEqual(dict(d.properties)["Blocks orthogonal to the quadratic model"], "no")
        self.assertFalse(block_dummies_orthogonal(d))

    def test_fractional_cube(self):
        d = doe.build_rsm(numeric(6), "ccd", randomize=False)
        self.assertEqual((d.fraction.runs, d.fraction.resolution), (32, 6))
        d = doe.build_rsm(numeric(6), "ccd", fraction=False, randomize=False)
        self.assertEqual(d.fraction.runs, 64)
        self.assertEqual(doe.ccd_cube(4).runs, 16)          # full below 5 factors

    def test_errors(self):
        for kwargs, msg in (({"alpha": "orthogonal"}, "needs blocks"),
                            ({"alpha": "custom"}, "positive"),
                            ({"blocks": 4}, "replicates times"),
                            ({"blocks": 5}, "confound two-factor")):
            with self.assertRaisesRegex(doe.DesignError, msg):
                doe.build_rsm(numeric(4), "ccd", **kwargs)
        with self.assertRaisesRegex(doe.DesignError, "2 to 7"):
            doe.build_rsm(numeric(8), "ccd")
        with self.assertRaisesRegex(doe.DesignError, "text levels"):
            doe.build_rsm(numeric(2) + [doe.Factor("Cat", "", "a", "b")], "ccd")
        d = doe.build_rsm(numeric(2), "ccd", centre_points=0, randomize=False)
        self.assertTrue(any("lack of fit" in n for n in d.notes))


def block_dummies_orthogonal(d):
    """Independent check: centred block indicators are orthogonal to the model columns."""
    X = np.array([quad(r.coded)[1:] for r in d.runs])
    blocks = sorted({r.block for r in d.runs})
    Z = np.array([[1.0 if r.block == b else 0.0 for b in blocks] for r in d.runs])
    Z -= Z.mean(axis=0)
    return np.allclose(Z.T @ X, 0, atol=1e-9)


class TestBoxBehnken(unittest.TestCase):
    def test_runs_levels_and_estimability(self):
        for k, runs in ((3, 15), (4, 27), (5, 46), (6, 54), (7, 62)):   # Box & Behnken (1960)
            d = doe.build_rsm(numeric(k), "bbd", randomize=False)
            self.assertEqual(len(d.runs), runs)
            x = np.array([r.coded for r in d.runs])
            self.assertTrue(all(set(x[:, i]) == {-1, 0, 1} for i in range(k)))
            self.assertTrue(np.all(x.sum(axis=0) == 0))
            X = np.array([quad(p) for p in x])
            self.assertEqual(np.linalg.matrix_rank(X), doe.model_terms(k))
            self.assertFalse(np.any(np.all(np.abs(x) == 1, axis=1)))   # no corner points

    def test_rotatable_only_for_4_and_7(self):
        for k in range(3, 8):
            d = doe.build_rsm(numeric(k), "bbd", randomize=False)
            self.assertEqual(dict(d.properties)["Rotatable"], "yes" if k in (4, 7) else "no")

    def test_seven_factor_sets_are_balanced(self):
        pairs = Counter(p for s in doe.BB_SETS[7] for p in combinations(sorted(s), 2))
        self.assertEqual(set(pairs.values()), {1})
        self.assertEqual(len(pairs), 21)

    def test_four_factors_in_three_orthogonal_blocks(self):
        d = doe.build_rsm(numeric(4), "bbd", blocks=3, randomize=False)
        self.assertEqual(len(d.runs), 3 * (8 + 2))
        self.assertEqual(dict(d.properties)["Blocks orthogonal to the quadratic model"], "yes")
        self.assertTrue(block_dummies_orthogonal(d))
        with self.assertRaisesRegex(doe.DesignError, "replicates times 1"):
            doe.build_rsm(numeric(5), "bbd", blocks=3)


class TestReadFactors(unittest.TestCase):
    def test_layouts(self):
        fs = doe.read_factors([("Temp", "°C", 40, 80), ("Time", "", 10.0, 30.0),
                               ("", "", "", ""), ("Cat", "", "Pd", "Pt")])
        self.assertEqual([f.label() for f in fs], ["Temp (°C)", "Time", "Cat"])
        self.assertFalse(fs[2].numeric)
        fs = doe.read_factors([("A", 1, 2), ("B", 3, 4)])
        self.assertEqual((fs[1].unit, fs[1].low), ("", 3))

    def test_errors(self):
        for rows, msg in (([("A", "", 1, 1), ("B", "", 0, 1)], "equal"),
                          ([("A", "", 1, ""), ("B", "", 0, 1)], "both a low"),
                          ([("A", "", 1, "x"), ("B", "", 0, 1)], "both be numbers"),
                          ([("A", "", 0, 1), ("A", "", 0, 1)], "different"),
                          ([("A", "", 0, 1)], "at least 2"),
                          ([("A", 0), ("B", 1)], "columns")):
            with self.assertRaisesRegex(doe.DesignError, msg):
                doe.read_factors(rows)


if __name__ == "__main__":
    unittest.main()
