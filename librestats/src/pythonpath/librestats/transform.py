"""Transformations, scaling, ranks and normal scores (pure Python).

Each function takes a list of floats with None for missing values and returns a
Transformed result whose `values` keep None in the same positions.

References
----------
Box, G. E. P., & Cox, D. R. (1964). An analysis of transformations. JRSS B, 26(2), 211-252.
Yeo, I.-K., & Johnson, R. A. (2000). A new family of power transformations to
    improve normality or symmetry. Biometrika, 87(4), 954-959.
Smithson, M., & Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood
    regression with beta-distributed dependent variables. Psychol. Methods, 11(1), 54-71.
van den Berg, R. A., Hoefsloot, H. C. J., Westerhuis, J. A., Smilde, A. K., &
    van der Werf, M. J. (2006). Centering, scaling, and transformations: improving
    the biological information content of metabolomics data. BMC Genomics, 7, 142.
Blom, G. (1958). Statistical Estimates and Transformed Beta-Variables. Wiley.
"""

import math
import statistics
from dataclasses import dataclass, field

from . import distributions


class TransformError(ValueError):
    """The transformation cannot be applied to these values."""


@dataclass
class Transformed:
    values: list
    params: list = field(default_factory=list)   # [(label, value), ...] for the report
    notes: list = field(default_factory=list)


def _present(values):
    return [v for v in values if v is not None]


def _apply(values, f):
    return [None if v is None else f(v) for v in values]


def _require(xs, minimum, what):
    if len(xs) < minimum:
        raise TransformError(f"{what} needs at least {minimum} values (got {len(xs)}).")


# --- Offsets ----------------------------------------------------------------------

OFFSETS = ("none", "one", "half_min")


def offset_value(values, rule):
    """Constant added before log / square root.

    "none" -> 0, "one" -> 1, "half_min" -> half the smallest positive value
    (common for zeros in concentration data), or a number.
    """
    if rule == "none":
        return 0.0
    if rule == "one":
        return 1.0
    if rule == "half_min":
        positive = [v for v in _present(values) if v > 0]
        if not positive:
            raise TransformError("There are no positive values to base the offset on.")
        return min(positive) / 2.0
    return float(rule)


def _offset_params(c):
    return [("Offset added", c)] if c else []


# --- log, square root, logit ------------------------------------------------------------

LOG_BASES = {"e": math.log, "10": math.log10, "2": math.log2}


def log_transform(values, base="e", offset="none"):
    c = offset_value(values, offset)
    bad = sum(1 for v in _present(values) if v + c <= 0)
    if bad:
        raise TransformError(f"{bad} value(s) are zero or negative after adding the offset "
                             f"({c:g}); choose a (larger) offset.")
    log = LOG_BASES[base]
    return Transformed(_apply(values, lambda v: log(v + c)),
                       [("Base", "e (ln)" if base == "e" else base)] + _offset_params(c))


def sqrt_transform(values, offset="none"):
    c = offset_value(values, offset)
    bad = sum(1 for v in _present(values) if v + c < 0)
    if bad:
        raise TransformError(f"{bad} value(s) are negative after adding the offset ({c:g}).")
    return Transformed(_apply(values, lambda v: math.sqrt(v + c)), _offset_params(c))


def logit_transform(values, percent=False, squeeze=True):
    """ln(p / (1 - p)) for proportions; percentages are divided by 100 first.

    If any value is exactly 0 or 1 and `squeeze` is set, all values are
    compressed with p' = (p (n - 1) + 0.5) / n (Smithson & Verkuilen, 2006).
    """
    scale = 100.0 if percent else 1.0
    props = _apply(values, lambda v: v / scale)
    xs = _present(props)
    _require(xs, 1, "Logit")
    outside = sum(1 for p in xs if p < 0 or p > 1)
    if outside:
        unit = "0-100" if percent else "0-1"
        raise TransformError(f"{outside} value(s) are outside the range {unit}.")
    params = [("Input", "percentages (divided by 100)" if percent else "proportions")]
    notes = []
    if any(p in (0.0, 1.0) for p in xs):
        if not squeeze:
            raise TransformError("Some values are exactly 0 or 1, where the logit is infinite. "
                                 "Allow squeezing, or remove these values.")
        n = len(xs)
        props = _apply(props, lambda p: (p * (n - 1) + 0.5) / n)
        params.append(("Squeeze n", n))
        notes.append("Values of 0 or 1 present: all values squeezed with (p(n-1) + 0.5)/n "
                     "(Smithson & Verkuilen, 2006).")
    return Transformed(_apply(props, lambda p: math.log(p / (1.0 - p))), params, notes)


# --- One-parameter maximisation ---------------------------------------------------------

_GOLDEN = (math.sqrt(5) - 1) / 2


def _maximise(f, lo, hi, tol=1e-10):
    """Golden-section search for the maximum of a unimodal f on [lo, hi]."""
    a, b = lo, hi
    c, d = b - _GOLDEN * (b - a), a + _GOLDEN * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - _GOLDEN * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + _GOLDEN * (b - a)
            fd = f(d)
    return 0.5 * (a + b)


def _bisect(f, a, b, tol=1e-12):
    """Root of f between a and b (f(a), f(b) of opposite sign)."""
    fa = f(a)
    for _ in range(200):
        m = 0.5 * (a + b)
        if abs(b - a) < tol:
            break
        fm = f(m)
        if (fm < 0) == (fa < 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


LAMBDA_RANGE = (-10.0, 10.0)


def _profile(llf, alpha):
    """ML lambda and likelihood-ratio confidence interval."""
    lo, hi = LAMBDA_RANGE
    lam = _maximise(llf, lo, hi)
    cut = llf(lam) - 0.5 * distributions.chi2_isf(alpha, 1)
    g = lambda t: llf(t) - cut
    ci_lo = _bisect(g, lo, lam) if g(lo) < 0 else None
    ci_hi = _bisect(g, lam, hi) if g(hi) < 0 else None
    notes = []
    if min(lam - lo, hi - lam) < 1e-6:
        notes.append(f"Lambda reached the search limit ({lo:g} to {hi:g}); "
                     "the data may not suit this transformation.")
    return lam, ci_lo, ci_hi, notes


def _var(ys):
    return statistics.pvariance(ys)


# --- Box-Cox ------------------------------------------------------------------------

def _boxcox(x, lam):
    if abs(lam) < 1e-12:
        return math.log(x)
    return math.expm1(lam * math.log(x)) / lam


def boxcox_llf(lam, xs):
    """Profile log-likelihood of Box-Cox lambda (same as scipy.stats.boxcox_llf)."""
    ys = [_boxcox(x, lam) for x in xs]
    var = _var(ys)
    if var <= 0:
        return -math.inf
    return (lam - 1) * sum(math.log(x) for x in xs) - len(xs) / 2 * math.log(var)


def boxcox(values, lmbda=None, offset="none", alpha=0.05):
    """Box-Cox transform; lambda estimated by maximum likelihood unless given."""
    c = offset_value(values, offset)
    xs = [v + c for v in _present(values)]
    _require(xs, 3, "Box-Cox")
    if min(xs) <= 0:
        raise TransformError("Box-Cox needs positive values; choose an offset or use "
                             "Yeo-Johnson.")
    if len(set(xs)) < 2:
        raise TransformError("All values are equal.")
    params, notes = _offset_params(c), []
    if lmbda is None:
        lam, lo, hi, notes = _profile(lambda t: boxcox_llf(t, xs), alpha)
        params += [("Lambda (ML)", lam),
                   (f"{100 * (1 - alpha):g}% CI lower", lo if lo is not None else "below range"),
                   (f"{100 * (1 - alpha):g}% CI upper", hi if hi is not None else "above range")]
    else:
        lam = float(lmbda)
        params.append(("Lambda (fixed)", lam))
    return Transformed(_apply(values, lambda v: _boxcox(v + c, lam)), params, notes)


# --- Yeo-Johnson ----------------------------------------------------------------------

def _yeojohnson(x, lam):
    if x >= 0:
        if abs(lam) < 1e-12:
            return math.log1p(x)
        return math.expm1(lam * math.log1p(x)) / lam
    if abs(lam - 2) < 1e-12:
        return -math.log1p(-x)
    return -math.expm1((2 - lam) * math.log1p(-x)) / (2 - lam)


def yeojohnson_llf(lam, xs):
    """Profile log-likelihood of Yeo-Johnson lambda (same as scipy.stats.yeojohnson_llf)."""
    ys = [_yeojohnson(x, lam) for x in xs]
    var = _var(ys)
    if var <= 0:
        return -math.inf
    jac = sum(math.copysign(math.log1p(abs(x)), x) for x in xs)
    return -len(xs) / 2 * math.log(var) + (lam - 1) * jac


def yeojohnson(values, lmbda=None, alpha=0.05):
    """Yeo-Johnson transform (any real values); lambda by maximum likelihood unless given."""
    xs = _present(values)
    _require(xs, 3, "Yeo-Johnson")
    if len(set(xs)) < 2:
        raise TransformError("All values are equal.")
    params, notes = [], []
    if lmbda is None:
        lam, lo, hi, notes = _profile(lambda t: yeojohnson_llf(t, xs), alpha)
        params += [("Lambda (ML)", lam),
                   (f"{100 * (1 - alpha):g}% CI lower", lo if lo is not None else "below range"),
                   (f"{100 * (1 - alpha):g}% CI upper", hi if hi is not None else "above range")]
    else:
        lam = float(lmbda)
        params.append(("Lambda (fixed)", lam))
    return Transformed(_apply(values, lambda v: _yeojohnson(v, lam)), params, notes)


# --- Scaling (van den Berg et al., 2006) ------------------------------------------------

SCALINGS = ("centre", "z", "range", "pareto")


def scale(values, method):
    """centre: x - mean; z (autoscaling): (x - mean) / sd; range: (x - min) / (max - min);
    pareto: (x - mean) / sqrt(sd). sd uses n - 1."""
    xs = _present(values)
    _require(xs, 2, "Scaling")
    mean = statistics.fmean(xs)
    if method == "centre":
        return Transformed(_apply(values, lambda v: v - mean), [("Mean", mean)])
    if method == "range":
        lo, hi = min(xs), max(xs)
        if hi == lo:
            raise TransformError("All values are equal (range is 0).")
        return Transformed(_apply(values, lambda v: (v - lo) / (hi - lo)),
                           [("Minimum", lo), ("Maximum", hi)])
    sd = statistics.stdev(xs)
    if sd == 0:
        raise TransformError("All values are equal (standard deviation is 0).")
    if method == "z":
        return Transformed(_apply(values, lambda v: (v - mean) / sd), [("Mean", mean), ("SD", sd)])
    if method == "pareto":
        root = math.sqrt(sd)
        return Transformed(_apply(values, lambda v: (v - mean) / root),
                           [("Mean", mean), ("SD", sd), ("sqrt(SD)", root)])
    raise ValueError(f"Unknown scaling {method!r}")


# --- Ranks and normal scores ------------------------------------------------------------

TIES = ("average", "min", "max", "dense")


def ranks(values, ties="average"):
    """Ranks (1 = smallest) with ties handled like scipy.stats.rankdata."""
    order = sorted((v, i) for i, v in enumerate(values) if v is not None)
    out = [None] * len(values)
    dense = 0
    k = 0
    while k < len(order):
        j = k
        while j + 1 < len(order) and order[j + 1][0] == order[k][0]:
            j += 1
        dense += 1
        rank = {"average": (k + j) / 2 + 1, "min": k + 1, "max": j + 1, "dense": dense}[ties]
        for t in range(k, j + 1):
            out[order[t][1]] = float(rank)
        k = j + 1
    return Transformed(out, [("Ties", ties)])


# a in Phi^-1((r - a) / (n + 1 - 2a))
NORMAL_SCORES = {"blom": 3 / 8, "vdw": 0.0, "tukey": 1 / 3, "rankit": 0.5}
NORMAL_SCORE_NAMES = {"blom": "Blom", "vdw": "van der Waerden", "tukey": "Tukey",
                      "rankit": "Rankit"}


def normal_scores(values, method="blom"):
    """Rank-based inverse normal scores; ties get the score of their average rank."""
    a = NORMAL_SCORES[method]
    r = ranks(values).values
    n = sum(1 for v in r if v is not None)
    _require([v for v in r if v is not None], 1, "Normal scores")
    out = _apply(r, lambda rank: distributions.norm_ppf((rank - a) / (n + 1 - 2 * a)))
    return Transformed(out, [("Method", f"{NORMAL_SCORE_NAMES[method]}: "
                                        f"Phi^-1((r - {a:.4g}) / (n + {1 - 2 * a:.4g}))"),
                             ("n", n)])
