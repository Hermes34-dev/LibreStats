"""Outlier detection methods (pure Python).

Univariate methods take a list of numbers and return a Result whose `flags`
point back to positions in that list. Mahalanobis takes rows of numbers.

References
----------
Tukey, J. W. (1977). Exploratory Data Analysis. Addison-Wesley.
Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting
    outliers: Do not use standard deviation around the mean, use absolute
    deviation around the median. J. Exp. Soc. Psychol., 49(4), 764-766.
Grubbs, F. E. (1969). Procedures for detecting outlying observations in
    samples. Technometrics, 11(1), 1-21.
Dixon, W. J. (1951). Ratios involving extreme values. Ann. Math. Stat.,
    22(1), 68-78.
Mahalanobis, P. C. (1936). On the generalised distance in statistics.
    Proc. Natl. Inst. Sci. India, 2, 49-55.
"""

import math
import statistics
from dataclasses import dataclass, field

from . import distributions, linalg
from .dixon_table import ALPHAS as DIXON_ALPHAS, CRITICAL as DIXON_CRITICAL

# Consistency constant for MAD under normality: 1 / Phi^-1(3/4) (Leys et al., 2013).
MAD_CONSTANT = 1.4826


class OutlierError(ValueError):
    """The method cannot be applied to this data (too few values, zero spread, ...)."""


@dataclass
class Flag:
    index: int      # position in the input list (row position for Mahalanobis)
    value: object   # the flagged value (a tuple of values for Mahalanobis)
    score: float    # method-specific statistic for this observation
    side: str = ""  # "low" / "high" for univariate methods


@dataclass
class Result:
    method: str
    n: int
    summary: list = field(default_factory=list)   # [(label, value), ...] for the report
    flags: list = field(default_factory=list)     # [Flag, ...]
    steps: list = field(default_factory=list)     # per-step rows (Grubbs iterations, Dixon)
    scores: list = field(default_factory=list)    # per-observation scores (Mahalanobis: (D², p))
    notes: list = field(default_factory=list)


def _require(values, minimum, method):
    if len(values) < minimum:
        raise OutlierError(f"{method} needs at least {minimum} values (got {len(values)}).")


def quantile(sorted_values, p):
    """Quantile by linear interpolation (Hyndman & Fan type 7).

    Same as Calc's QUARTILE / PERCENTILE (and QUARTILE.INC / PERCENTILE.INC).
    """
    n = len(sorted_values)
    h = (n - 1) * p
    lo = math.floor(h)
    hi = min(lo + 1, n - 1)
    return sorted_values[lo] + (h - lo) * (sorted_values[hi] - sorted_values[lo])


# --- Tukey fences -------------------------------------------------------------

def tukey_fences(values, k=1.5):
    """Flag values outside [Q1 - k*IQR, Q3 + k*IQR].

    Score = distance beyond the nearer quartile, in IQR units.
    """
    _require(values, 3, "Tukey fences")
    s = sorted(values)
    q1, q3 = quantile(s, 0.25), quantile(s, 0.75)
    iqr = q3 - q1
    lower, upper = q1 - k * iqr, q3 + k * iqr
    res = Result("Tukey fences", len(values), summary=[
        ("Q1", q1), ("Q3", q3), ("IQR", iqr), ("k", k),
        ("Lower fence", lower), ("Upper fence", upper)])
    for i, x in enumerate(values):
        if x < lower or x > upper:
            side = "low" if x < lower else "high"
            dist = (q1 - x) if side == "low" else (x - q3)
            res.flags.append(Flag(i, x, dist / iqr if iqr else math.inf, side))
    if iqr == 0:
        res.notes.append("IQR is 0: every value different from the quartiles is flagged.")
    return res


# --- Robust z-scores (median / MAD) -------------------------------------------

def mad_zscores(values, threshold=2.5):
    """Flag values with |x - median| / MAD > threshold, MAD = 1.4826 * median|x - median|.

    Leys et al. (2013) suggest 2.5 (moderately conservative); 3 is very
    conservative and 2 poorly conservative.
    """
    _require(values, 3, "Robust z-scores (MAD)")
    med = statistics.median(values)
    raw_mad = statistics.median(abs(x - med) for x in values)
    if raw_mad == 0:
        raise OutlierError("MAD is 0 (at least half of the values are identical), "
                           "so robust z-scores cannot be computed.")
    mad = MAD_CONSTANT * raw_mad
    res = Result("Robust z-scores (MAD)", len(values), summary=[
        ("Median", med), ("MAD (scaled, b = 1.4826)", mad), ("Threshold |z|", threshold),
        ("Lower limit", med - threshold * mad), ("Upper limit", med + threshold * mad)])
    for i, x in enumerate(values):
        z = (x - med) / mad
        if abs(z) > threshold:
            res.flags.append(Flag(i, x, z, "low" if z < 0 else "high"))
    return res


# --- Grubbs test --------------------------------------------------------------

def grubbs_critical(n, alpha):
    """Two-sided critical value of Grubbs' G for sample size n."""
    t = distributions.t_isf(alpha / (2 * n), n - 2)
    return (n - 1) / math.sqrt(n) * math.sqrt(t * t / (n - 2 + t * t))


def grubbs_pvalue(g, n):
    """Two-sided p-value for Grubbs' G (Bonferroni bound, exact enough for p < 0.5)."""
    denom = (n - 1) ** 2 - n * g * g
    if denom <= 0:
        return 0.0
    t = math.sqrt(n * (n - 2) * g * g / denom)
    return min(1.0, 2 * n * distributions.t_sf(t, n - 2))


def grubbs(values, alpha=0.05, iterate=False):
    """Two-sided Grubbs test for a single outlier; optionally repeated.

    Assumes the data (without outliers) are approximately normal. When
    `iterate` is True the flagged value is removed and the test repeated until
    no outlier is found or fewer than 3 values remain.
    """
    _require(values, 3, "Grubbs test")
    remaining = list(enumerate(values))
    res = Result("Grubbs test", len(values), summary=[("alpha (two-sided)", alpha)])
    while len(remaining) >= 3:
        xs = [x for _, x in remaining]
        n = len(xs)
        mean = statistics.fmean(xs)
        sd = statistics.stdev(xs)
        if sd == 0:
            res.notes.append("Standard deviation is 0: no further testing possible.")
            break
        pos = max(range(n), key=lambda j: abs(xs[j] - mean))
        index, x = remaining[pos]
        g = abs(x - mean) / sd
        g_crit = grubbs_critical(n, alpha)
        p = grubbs_pvalue(g, n)
        outlier = g > g_crit
        res.steps.append({"n": n, "mean": mean, "sd": sd, "index": index, "value": x,
                          "G": g, "G_crit": g_crit, "p": p, "outlier": outlier})
        if not outlier:
            break
        res.flags.append(Flag(index, x, g, "low" if x < mean else "high"))
        if not iterate:
            break
        del remaining[pos]
    if iterate:
        res.notes.append("Repeated Grubbs testing can miss outliers that mask each other.")
    return res


# --- Dixon test ---------------------------------------------------------------

# rIJ as (gap, trim): numerator spans `gap` order statistics, denominator skips
# `trim` values at the opposite end. Must match tools/make_dixon_table.py.
_DIXON_SPEC = {"r10": (1, 0), "r11": (1, 1), "r21": (2, 1), "r22": (2, 2)}


def _ratio(num, den):
    return num / den if den > 0 else 0.0


def dixon(values, alpha=0.05):
    """Two-sided Dixon test of the most extreme value (n = 3..30)."""
    n = len(values)
    if not 3 <= n <= 30:
        raise OutlierError(f"Dixon's test needs 3 to 30 values (got {n}).")
    if alpha not in DIXON_ALPHAS:
        raise OutlierError(f"Dixon's test supports alpha = {DIXON_ALPHAS}.")
    kind, crits = DIXON_CRITICAL[n]
    gap, trim = _DIXON_SPEC[kind]
    order = sorted(range(n), key=lambda i: values[i])
    x = [values[i] for i in order]
    if x[n - 1 - trim] == x[0] and x[n - 1] == x[trim]:
        raise OutlierError("Too many tied values: Dixon's ratio cannot be computed.")
    r_low = _ratio(x[gap] - x[0], x[n - 1 - trim] - x[0])
    r_high = _ratio(x[n - 1] - x[n - 1 - gap], x[n - 1] - x[trim])
    side = "high" if r_high >= r_low else "low"
    r = max(r_low, r_high)
    index = order[-1] if side == "high" else order[0]
    crit = crits[DIXON_ALPHAS.index(alpha)]

    # p-value bracket from the tabulated levels
    below = [a for a, c in zip(DIXON_ALPHAS, crits) if r > c]
    p_text = f"p < {min(below):g}" if below else f"p > {max(DIXON_ALPHAS):g}"

    res = Result("Dixon test", n, summary=[("alpha (two-sided)", alpha), ("Ratio", kind)])
    res.steps.append({"n": n, "ratio": kind, "r_low": r_low, "r_high": r_high,
                      "index": index, "value": values[index], "r": r, "r_crit": crit,
                      "p": p_text, "outlier": r > crit})
    if r > crit:
        res.flags.append(Flag(index, values[index], r, side))
    return res


# --- Mahalanobis distance -----------------------------------------------------

def mahalanobis(rows, alpha=0.001):
    """Flag rows whose squared Mahalanobis distance exceeds the chi-square cutoff.

    D^2 = (x - mean)^T S^-1 (x - mean), with the sample mean and covariance
    (n - 1 denominator) of all rows. Under multivariate normality D^2 is
    approximately chi-square with p degrees of freedom; alpha = 0.001 is the
    usual convention (Tabachnick & Fidell, 2013). Rows must be complete.
    """
    if not rows:
        raise OutlierError("Mahalanobis distance needs data.")
    p = len(rows[0])
    n = len(rows)
    if p < 2:
        raise OutlierError("Mahalanobis distance needs at least 2 variables.")
    if n < p + 2:
        raise OutlierError(f"Mahalanobis distance with {p} variables needs at least "
                           f"{p + 2} complete rows (got {n}).")
    means = [statistics.fmean(r[j] for r in rows) for j in range(p)]
    centred = [[r[j] - means[j] for j in range(p)] for r in rows]
    cov = [[sum(c[i] * c[j] for c in centred) / (n - 1) for j in range(p)] for i in range(p)]
    try:
        L = linalg.cholesky(cov)
    except linalg.SingularMatrixError:
        raise OutlierError("The covariance matrix is singular: a variable is constant "
                           "or a linear combination of the others.") from None
    cutoff = distributions.chi2_isf(alpha, p)
    res = Result("Mahalanobis distance", n, summary=[
        ("Variables (df)", p), ("alpha", alpha), ("Critical D²", cutoff)])
    max_d2 = (n - 1) ** 2 / n  # largest D² any row of this sample can have
    if cutoff >= max_d2:
        res.notes.append(f"With n = {n}, no row can exceed the critical value (maximum "
                         f"possible D² is {max_d2:.3g}); use more data or a larger alpha.")
    for i, c in enumerate(centred):
        d2 = sum(v * v for v in linalg.forward_solve(L, c))
        pval = distributions.chi2_sf(d2, p)
        res.scores.append((d2, pval))
        if d2 > cutoff:
            res.flags.append(Flag(i, tuple(rows[i]), d2))
    return res
