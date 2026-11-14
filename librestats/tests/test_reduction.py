"""Check librestats.reduction against independent implementations with numpy and
statsmodels: backward elimination by p-value and all-subsets search by AICc, with
and without model hierarchy, locked terms, manual removal and the reduced cubic."""

import math
import os
import sys
import unittest
from itertools import combinations, combinations_with_replacement

import numpy as np
import statsmodels.api as sm
from scipy.stats import f as f_dist

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import doe, reduction, rsm  # noqa: E402

LETTERS = "ABCDEFGH"


def factors(k):
    return [doe.Factor(f"X{i}", "", 0.0, 1.0) for i in range(k)]


def design_data(kind, k, seed=0, model_y=None):
    d = doe.build_rsm(factors(k), kind, randomize=False)
    x = [r.coded for r in d.runs]
    rng = np.random.default_rng(seed)
    f = model_y or (lambda v: 10 + 2 * v[0] - v[1] ** 2 + 0.5 * v[0] * v[1])
    y = [float(f(v) + rng.normal()) for v in x]
    return x, [r.block for r in d.runs], y


# --- independent reference implementation ---------------------------------------------

def monos(k, model):
    out = [(i,) for i in range(k)]
    if model in ("quadratic",):
        out += [(i, i) for i in range(k)] + list(combinations(range(k), 2))
    return out


def name(m):
    s = ""
    for i in sorted(set(m)):
        c = m.count(i)
        s += LETTERS[i] + (f"^{c}" if c > 1 else "")
    return s


def ref_parents(m):
    return {tuple(sorted(sub)) for size in range(1, len(m)) for sub in combinations(m, size)}


def design_matrix(x, terms):
    return np.array([[1.0] + [np.prod([v[i] for i in m]) for m in terms] for v in x])


def ref_aicc(x, y, terms):
    X = design_matrix(x, terms)
    n, p = X.shape
    if n - p - 1 <= 0:
        return None
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    sse = float(((y - X @ beta) ** 2).sum())
    return n * math.log(sse / n) + 2 * p + 2 * p * (p + 1) / (n - p - 1)


def ref_backward(x, y, terms, alpha, hierarchy=True, locked=()):
    current, steps = list(terms), []
    while True:
        res = sm.OLS(np.array(y), design_matrix(x, current)).fit()
        pv = dict(zip(current, res.pvalues[1:]))
        options = [m for m in current if m not in locked and
                   (not hierarchy or not any(m in ref_parents(u) for u in current))]
        if not options:
            return steps, current
        worst = max(options, key=lambda m: pv[m])
        if pv[worst] <= alpha:
            return steps, current
        steps.append((name(worst), pv[worst]))
        current.remove(worst)


def ref_subsets(x, y, terms, hierarchy=True):
    out = []
    for r in range(len(terms) + 1):
        for sub in combinations(terms, r):
            if hierarchy and any(not (ref_parents(m) <= set(sub)) for m in sub):
                continue
            a = ref_aicc(x, np.array(y), list(sub))
            if a is not None:
                out.append((a, len(sub), sorted(name(m) for m in sub)))
    out.sort(key=lambda t: (t[0], t[1]))
    return out


# --- tests ---------------------------------------------------------------------------

class TestHierarchy(unittest.TestCase):
    def test_parents(self):
        self.assertEqual(reduction.parents("A^2B"), {"A", "B", "A^2", "AB"})
        self.assertEqual(reduction.parents("ABC"), {"A", "B", "C", "AB", "AC", "BC"})
        self.assertEqual(reduction.parents("A^3"), {"A", "A^2"})
        self.assertEqual(reduction.parents("B"), set())

    def test_removable_and_missing_parents(self):
        start = ["A", "B", "A^2", "B^2", "AB"]
        self.assertEqual(reduction.removable(start, start, set()), ["A^2", "B^2", "AB"])
        self.assertEqual(reduction.removable(["A", "B", "B^2"], start, {"B^2"}), ["A"])
        self.assertEqual(reduction.removable(start, start, set(), hierarchy=False),
                         ["A", "B", "A^2", "B^2", "AB"])
        self.assertEqual(reduction.missing_parents(["A", "AB"], start), {"AB": ["B"]})
        # in a reduced cubic only parents present in the starting model count
        self.assertEqual(reduction.missing_parents(["A", "B", "A^2", "A^2B"],
                                                   ["A", "B", "A^2", "A^2B"]),
                         {})


class TestBackward(unittest.TestCase):
    def test_against_statsmodels(self):
        for kind, k, seed in (("ccd", 3, 0), ("ccd", 3, 5), ("bbd", 4, 1), ("ccd", 4, 2)):
            x, blocks, y = design_data(kind, k, seed)
            for alpha in (0.05, 0.10):
                for hierarchy in (True, False):
                    red = reduction.reduce(x, blocks, y, "quadratic", "backward", alpha,
                                           hierarchy)
                    steps, final = ref_backward(x, y, monos(k, "quadratic"), alpha,
                                                hierarchy)
                    msg = (kind, k, seed, alpha, hierarchy)
                    self.assertEqual([t for t, _ in red.steps], [t for t, _ in steps], msg)
                    np.testing.assert_allclose([p for _, p in red.steps],
                                               [p for _, p in steps], rtol=1e-8)
                    self.assertEqual(sorted(t.name for t in red.final.terms[1:]),
                                     sorted(name(m) for m in final), msg)
                    if hierarchy:
                        self.assertEqual(reduction.missing_parents(
                            [t.name for t in red.final.terms[1:]],
                            [name(m) for m in monos(k, "quadratic")]), {})

    def test_locked_terms_stay(self):
        x, blocks, y = design_data("ccd", 3, 0)
        red = reduction.reduce(x, blocks, y, "quadratic", "backward", 0.05,
                               locked=["C^2"])
        names = [t.name for t in red.final.terms]
        self.assertIn("C^2", names)
        self.assertIn("C", names)                    # parent of a locked term
        steps, final = ref_backward(x, y, monos(3, "quadratic"), 0.05,
                                    locked=[(2, 2)])
        self.assertEqual([t for t, _ in red.steps], [t for t, _ in steps])

    def test_reduced_cubic_start(self):
        x, blocks, y = design_data("cccd", 3, 3)
        red = reduction.reduce(x, blocks, y, "cubic", "backward", 0.05)
        start = [t.name for t in red.start.terms[1:]]
        self.assertEqual(len(start), 16)             # rank 17 of 20, minus the constant
        self.assertEqual(reduction.missing_parents([t.name for t in red.final.terms[1:]],
                                                   start), {})
        for t, p in red.steps:
            self.assertGreater(p, 0.05)


class TestSubsets(unittest.TestCase):
    def test_against_brute_force(self):
        for kind, k, seed, hierarchy in (("ccd", 3, 0, True), ("ccd", 3, 4, False),
                                         ("bbd", 3, 2, True), ("ccd", 4, 1, True)):
            x, blocks, y = design_data(kind, k, seed)
            red = reduction.reduce(x, blocks, y, "quadratic", "subsets",
                                   hierarchy=hierarchy, top=5)
            ref = ref_subsets(x, y, monos(k, "quadratic"), hierarchy)
            msg = (kind, k, seed, hierarchy)
            self.assertEqual(red.evaluated, len(ref), msg)
            np.testing.assert_allclose([c.aicc for c in red.candidates],
                                       [a for a, _, _ in ref[:5]], rtol=1e-9, err_msg=str(msg))
            self.assertEqual([sorted(c.terms) for c in red.candidates],
                             [t for _, _, t in ref[:5]], msg)
            self.assertEqual(sorted(t.name for t in red.final.terms[1:]), ref[0][2], msg)

    def test_statistics_of_the_shown_models(self):
        x, blocks, y = design_data("ccd", 3, 0)
        red = reduction.reduce(x, blocks, y, "quadratic", "subsets", top=4, rank=2)
        yy = np.array(y)
        groups = {}
        for v, val in zip(x, y):
            groups.setdefault(tuple(v), []).append(val)
        pure_ss = sum(((np.array(g) - np.mean(g)) ** 2).sum() for g in groups.values())
        pure_df = sum(len(g) - 1 for g in groups.values())
        for c in red.candidates:
            terms = [rsm.name_monomial(t) for t in c.terms]
            res = sm.OLS(yy, design_matrix(x, terms)).fit()
            self.assertAlmostEqual(c.r2, res.rsquared, places=10)
            self.assertAlmostEqual(c.r2_adj, res.rsquared_adj, places=10)
            h = res.get_influence().hat_matrix_diag
            press = ((res.resid / (1 - h)) ** 2).sum()
            self.assertAlmostEqual(c.r2_pred, 1 - press / ((yy - yy.mean()) ** 2).sum(),
                                   places=10)
            lof_df = res.df_resid - pure_df
            f = ((res.ssr - pure_ss) / lof_df) / (pure_ss / pure_df)
            self.assertAlmostEqual(c.lof_p, f_dist.sf(f, lof_df, pure_df), places=10)
        self.assertEqual(red.rank, 2)
        self.assertEqual(sorted(t.name for t in red.final.terms[1:]),
                         sorted(red.candidates[1].terms))

    def test_guards(self):
        x, blocks, y = design_data("ccd", 6, 0)            # 27 candidate terms
        with self.assertRaisesRegex(reduction.ReductionError, "too many"):
            reduction.reduce(x, blocks, y, "quadratic", "subsets")
        x, blocks, y = design_data("ccd", 3, 0)
        with self.assertRaisesRegex(reduction.ReductionError, "rank 500 is not available"):
            reduction.reduce(x, blocks, y, "quadratic", "subsets", top=5, rank=500)
        red = reduction.reduce(x, blocks, y, "quadratic", "subsets", top=3, rank=7)
        self.assertEqual(len(red.candidates), 3)
        self.assertEqual(sorted(t.name for t in red.final.terms[1:]), sorted(red.chosen.terms))


class TestManualAndWarnings(unittest.TestCase):
    def test_manual(self):
        x, blocks, y = design_data("ccd", 3, 0)
        red = reduction.reduce(x, blocks, y, "quadratic", "manual", remove=["AC", "C^2"])
        self.assertEqual(red.removed, ["C^2", "AC"])
        self.assertNotIn("AC", [t.name for t in red.final.terms])
        with self.assertRaisesRegex(reduction.ReductionError,
                                    "C cannot be removed while C\\^2, AC and BC are in the model"):
            reduction.reduce(x, blocks, y, "quadratic", "manual", remove=["C"])
        with self.assertRaisesRegex(reduction.ReductionError, "locked"):
            reduction.reduce(x, blocks, y, "quadratic", "manual", remove=["AB"],
                             locked=["AB"])
        red = reduction.reduce(x, blocks, y, "quadratic", "manual", remove=["C"],
                               hierarchy=False)
        self.assertTrue(any("hierarchy is off" in w and "AC without C" in w
                            for w in red.warnings))

    def test_aliased_term_cannot_be_locked(self):
        x, blocks, y = design_data("ccd", 3, 0)
        with self.assertRaisesRegex(reduction.ReductionError, "aliased"):
            reduction.reduce(x, blocks, y, "cubic", "backward", locked=["A^3"])

    def test_summary_text(self):
        x, blocks, y = design_data("ccd", 3, 0)
        red = reduction.reduce(x, blocks, y, "quadratic", "backward", 0.05)
        text = red.summary("full quadratic")
        self.assertTrue(text.startswith("Reduced from full quadratic by backward elimination "
                                        "(alpha = 0.05; hierarchy kept). Removed: "), text)
        self.assertTrue(text.endswith(", ".join(t for t, _ in red.steps) + "."))


if __name__ == "__main__":
    unittest.main()


DESIGN2_Y = [24.28989323, 33.66123277, 39.83253039, 18.65134664, 12.52643918, 24.52400157,
             19.43370304, 18.04370825, 11.58758592, 19.25931338, 16.25878132, 9.714146726,
             13.70971019, 33.93868664, 22.13181921, 51.13647479, 26.76090827, 37.48920235,
             15.37816958, 31.05661974, 11.95990485, 18.31395824, 17.97089328, 11.73768574,
             8.318061785, 12.1663204, 14.31427388, 10.81211525]


def design2_without_run3():
    d = doe.build_rsm(factors(3), "cccd", randomize=False)
    x = [r.coded for i, r in enumerate(d.runs) if i != 2]
    y = [v for i, v in enumerate(DESIGN2_Y) if i != 2]
    return x, [1] * 27, y


class TestAliasMarksAndRegion(unittest.TestCase):
    def test_alias_marks(self):
        x, blocks, y = design2_without_run3()
        fit = rsm.fit_model(x, blocks, y, "cubic")
        self.assertEqual(rsm.alias_marks(fit.aliases),
                         {"AB^2": "aliased with AC^2", "A^2C": "aliased with B^2C",
                          "A^2B": "aliased with BC^2"})
        x2 = [r.coded for r in doe.build_rsm(factors(2), "ccd", randomize=False).runs]
        fit = rsm.fit_model(x2, [1] * len(x2), [float(i % 4) for i in range(len(x2))],
                            "cubic")
        # A^3 = 2 A - AB^2: AB^2 (same order) is marked, A (lower order) is not
        self.assertEqual(rsm.alias_marks(fit.aliases),
                         {"AB^2": "part of the A^3 alias", "A^2B": "part of the B^3 alias"})

    def test_corners_not_run(self):
        a = 2 ** 0.5
        ccd = [r.coded for r in doe.build_rsm(factors(2), "ccd", randomize=False).runs]
        self.assertEqual(len(rsm.corners_not_run([(-a, a), (-a, a)], ccd)), 4)
        cccd = [r.coded for r in doe.build_rsm(factors(3), "cccd", randomize=False).runs]
        b = 8 ** 0.25
        self.assertEqual(rsm.corners_not_run([(-b, b)] * 3, cccd), [])
        bbd = [r.coded for r in doe.build_rsm(factors(3), "bbd", randomize=False).runs]
        self.assertEqual(len(rsm.corners_not_run([(-1, 1)] * 3, bbd)), 8)

    def test_exclude_aliased(self):
        x, blocks, y = design2_without_run3()
        red = reduction.reduce(x, blocks, y, "cubic", "subsets", exclude_aliased=True)
        self.assertEqual(red.excluded, ["A^2B", "A^2C", "AB^2"])
        terms = {t.name for t in red.start.terms}
        self.assertFalse(terms & {"A^2B", "A^2C", "AB^2"})
        self.assertTrue(all(not (set(c.terms) & {"A^2B", "A^2C", "AB^2"})
                            for c in red.candidates))
        with self.assertRaisesRegex(reduction.ReductionError, "locked and excluded"):
            reduction.reduce(x, blocks, y, "cubic", "backward", locked=["AB^2"],
                             exclude_aliased=True)

    def test_users_case(self):
        # the user's run: 27 runs (run 3 removed), reduced cubic, all subsets
        x, blocks, y = design2_without_run3()
        red = reduction.reduce(x, blocks, y, "cubic", "subsets")
        self.assertEqual(red.start.full_terms, 20)
        self.assertEqual(red.start.full_terms - len(red.start.aliases), 17)
        self.assertEqual(len(red.final.terms), 9)
        self.assertEqual(sorted(red.chosen.terms),
                         sorted(["A", "B", "C", "B^2", "C^2", "AB", "BC", "AB^2"]))
        self.assertLess(red.final.minimum[1], min(y))      # predicted minimum < 8.32
