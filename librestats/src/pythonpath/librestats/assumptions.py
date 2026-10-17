"""Assumption checks: normality (Shapiro-Wilk) and equal variances (Levene).

Shapiro-Wilk follows Royston's algorithm AS R94 (Royston, 1995), which
approximates the coefficients and the null distribution of W for 3 <= n <= 5000
(Royston, 1992), as R's shapiro.test and scipy.stats.shapiro do. Expected normal
order statistics use the exact normal quantile, as in R.

Levene's test (Levene, 1960) is a one-way ANOVA on the absolute deviations
from each group's centre; centring on the median is the Brown-Forsythe
version (Brown & Forsythe, 1974), which is robust to non-normal data.
"""

import math
from dataclasses import dataclass, field

from .distributions import f_sf, norm_cdf, norm_ppf

SW_MAX_N = 5000


class AssumptionError(ValueError):
    pass


@dataclass
class TestResult:
    statistic: float
    p: float
    n: int
    params: list = field(default_factory=list)   # extra (label, value) pairs
    notes: list = field(default_factory=list)


def _poly(c, x):
    """c[0] + c[1] x + c[2] x^2 + ..."""
    return sum(ci * x ** i for i, ci in enumerate(c))


# --- Shapiro-Wilk ---------------------------------------------------------------------

_C1 = (0.0, 0.221157, -0.147981, -2.071190, 4.434685, -2.706056)
_C2 = (0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633)
_C3 = (0.5440, -0.39978, 0.025054, -6.714e-4)
_C4 = (1.3822, -0.77857, 0.062767, -0.0020322)
_C5 = (-1.5861, -0.31082, -0.083751, 0.0038915)
_C6 = (-0.4803, -0.082676, 0.0030302)
_G = (-2.273, 0.459)


def shapiro_coefficients(n):
    """Royston's approximation of the Shapiro-Wilk coefficients a_1..a_n (antisymmetric)."""
    if n == 3:
        upper = [math.sqrt(0.5)]
    else:
        half = n // 2
        # m_i for the upper half, largest first: m_n, m_(n-1), ...
        m = [-norm_ppf((i - 0.375) / (n + 0.25)) for i in range(1, half + 1)]
        summ2 = 2.0 * sum(v * v for v in m)
        ssumm2 = math.sqrt(summ2)
        rsn = 1.0 / math.sqrt(n)
        a1 = _poly(_C1, rsn) + m[0] / ssumm2
        if n > 5:
            a2 = m[1] / ssumm2 + _poly(_C2, rsn)
            fac = math.sqrt((summ2 - 2 * m[0] ** 2 - 2 * m[1] ** 2) /
                            (1 - 2 * a1 ** 2 - 2 * a2 ** 2))
            upper = [a1, a2] + [v / fac for v in m[2:]]
        else:
            fac = math.sqrt((summ2 - 2 * m[0] ** 2) / (1 - 2 * a1 ** 2))
            upper = [a1] + [v / fac for v in m[1:]]
    middle = [0.0] if n % 2 else []
    return [-v for v in upper] + middle + upper[::-1]


def shapiro_wilk(values):
    """Shapiro-Wilk W and its p-value (H0: the data come from a normal distribution)."""
    x = sorted(v for v in values if v is not None)
    n = len(x)
    if n < 3:
        raise AssumptionError("Shapiro-Wilk needs at least 3 values.")
    mean = sum(x) / n
    ss = sum((v - mean) ** 2 for v in x)
    if ss == 0 or x[0] == x[-1]:
        raise AssumptionError("All values are equal.")
    a = shapiro_coefficients(n)
    w = sum(ai * xi for ai, xi in zip(a, x)) ** 2 / ss
    w = min(w, 1.0)
    notes = []

    if n == 3:
        p = max(0.0, 6 / math.pi * (math.asin(math.sqrt(max(w, 0.75))) - math.pi / 3))
    else:
        w1 = math.log(1 - w) if w < 1 else -math.inf
        if n <= 11:
            gamma = _poly(_G, n)
            if w1 >= gamma:
                p = 0.0
            else:
                y = -math.log(gamma - w1)
                p = 1 - norm_cdf((y - _poly(_C3, n)) / math.exp(_poly(_C4, n)))
        else:
            ln = math.log(n)
            y = w1
            p = 1 - norm_cdf((y - _poly(_C5, ln)) / math.exp(_poly(_C6, ln)))
    if n > SW_MAX_N:
        notes.append(f"Royston's approximation is validated up to n = {SW_MAX_N}; "
                     "the p-value for larger samples is approximate.")
    if len(set(x)) < n:
        notes.append("Tied values are present; heavy rounding of the data can make the test "
                     "reject normality.")
    return TestResult(w, min(max(p, 0.0), 1.0), n, notes=notes)


# --- Levene / Brown-Forsythe ------------------------------------------------------------

def _median(xs):
    s = sorted(xs)
    k = len(s)
    return s[k // 2] if k % 2 else 0.5 * (s[k // 2 - 1] + s[k // 2])


def levene(groups, center="median"):
    """Levene's test of equal variances across groups (lists of numbers).
    center='median' is the Brown-Forsythe version, 'mean' the original."""
    if center not in ("median", "mean"):
        raise AssumptionError(f"Unknown centre {center!r}.")
    groups = [[v for v in g if v is not None] for g in groups]
    k = len(groups)
    if k < 2:
        raise AssumptionError("Levene's test needs at least 2 groups.")
    if any(len(g) < 2 for g in groups):
        raise AssumptionError("Every group needs at least 2 values.")
    centre = _median if center == "median" else (lambda g: sum(g) / len(g))
    z = [[abs(v - centre(g)) for v in g] for g in groups]
    n_total = sum(len(g) for g in groups)
    zbar = [sum(zi) / len(zi) for zi in z]
    zall = sum(sum(zi) for zi in z) / n_total
    between = sum(len(zi) * (zm - zall) ** 2 for zi, zm in zip(z, zbar))
    within = sum((v - zm) ** 2 for zi, zm in zip(z, zbar) for v in zi)
    df1, df2 = k - 1, n_total - k
    if within == 0:
        raise AssumptionError("The absolute deviations from the group centres do not vary "
                              "within any group, so the test cannot be computed.")
    f = (df2 / df1) * between / within
    return TestResult(f, f_sf(f, df1, df2), n_total, [("df1", df1), ("df2", df2)])
