"""Probability distributions in pure Python (Student's t, chi-square, F, normal).

Built on the regularized incomplete beta and gamma functions, evaluated with
continued fractions / series as in Press et al., Numerical Recipes (3rd ed.),
sections 6.2 and 6.4. Quantiles are found by bisection on the CDF.
"""

import math

_EPS = 1e-15
_TINY = 1e-300
_MAX_ITER = 1000


# --- Regularized incomplete beta I_x(a, b) ----------------------------------

def _beta_cf(a, b, x):
    """Continued fraction for I_x(a, b) (modified Lentz's method)."""
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > _TINY else _TINY)
    h = d
    for m in range(1, _MAX_ITER + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > _TINY else _TINY)
        c = 1.0 + aa / c
        c = c if abs(c) > _TINY else _TINY
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > _TINY else _TINY)
        c = 1.0 + aa / c
        c = c if abs(c) > _TINY else _TINY
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPS:
            return h
    raise ArithmeticError("incomplete beta continued fraction did not converge")


def betainc(a, b, x):
    """Regularized incomplete beta function I_x(a, b)."""
    if a <= 0 or b <= 0:
        raise ValueError("a and b must be positive")
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_front = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                 + a * math.log(x) + b * math.log1p(-x))
    front = math.exp(log_front)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_cf(a, b, x) / a
    return 1.0 - front * _beta_cf(b, a, 1.0 - x) / b


# --- Regularized incomplete gamma P(a, x), Q(a, x) ---------------------------

def _gamma_series(a, x):
    """P(a, x) by its series representation (good for x < a + 1)."""
    ap, total = a, 1.0 / a
    term = total
    for _ in range(_MAX_ITER):
        ap += 1.0
        term *= x / ap
        total += term
        if abs(term) < abs(total) * _EPS:
            return total * math.exp(-x + a * math.log(x) - math.lgamma(a))
    raise ArithmeticError("incomplete gamma series did not converge")


def _gamma_cf(a, x):
    """Q(a, x) by its continued fraction (good for x >= a + 1)."""
    b = x + 1.0 - a
    c = 1.0 / _TINY
    d = 1.0 / b
    h = d
    for i in range(1, _MAX_ITER + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = 1.0 / (d if abs(d) > _TINY else _TINY)
        c = b + an / c
        c = c if abs(c) > _TINY else _TINY
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPS:
            return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h
    raise ArithmeticError("incomplete gamma continued fraction did not converge")


def gammainc_upper(a, x):
    """Regularized upper incomplete gamma function Q(a, x) = 1 - P(a, x)."""
    if a <= 0:
        raise ValueError("a must be positive")
    if x <= 0.0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _gamma_series(a, x)
    return _gamma_cf(a, x)


# --- Quantile helper ----------------------------------------------------------

def _invert_decreasing(sf, p, lo=0.0):
    """Find x >= lo with sf(x) = p, for a decreasing survival function sf."""
    hi = max(1.0, lo + 1.0)
    while sf(hi) > p:
        lo, hi = hi, hi * 2.0
        if hi > 1e300:
            raise ArithmeticError("quantile out of range")
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if mid in (lo, hi):
            break
        if sf(mid) > p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --- Student's t --------------------------------------------------------------

def t_sf(t, df):
    """Upper tail probability P(T > t) for Student's t with df degrees of freedom."""
    if df <= 0:
        raise ValueError("df must be positive")
    tail = 0.5 * betainc(0.5 * df, 0.5, df / (df + t * t))
    return tail if t >= 0 else 1.0 - tail


def t_cdf(t, df):
    return t_sf(-t, df)


def t_isf(p, df):
    """Inverse survival function: t such that P(T > t) = p."""
    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    if p == 0.5:
        return 0.0
    if p > 0.5:
        return -t_isf(1.0 - p, df)
    return _invert_decreasing(lambda t: t_sf(t, df), p)


# --- Chi-square ---------------------------------------------------------------

def chi2_sf(x, df):
    """Upper tail probability P(X > x) for chi-square with df degrees of freedom."""
    if df <= 0:
        raise ValueError("df must be positive")
    return gammainc_upper(0.5 * df, 0.5 * x)


def chi2_cdf(x, df):
    return 1.0 - chi2_sf(x, df)


def chi2_isf(p, df):
    """Inverse survival function: x such that P(X > x) = p."""
    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    return _invert_decreasing(lambda x: chi2_sf(x, df), p)


# --- F ------------------------------------------------------------------------

def f_sf(x, df1, df2):
    """Upper tail probability P(F > x) for the F distribution with (df1, df2) df."""
    if df1 <= 0 or df2 <= 0:
        raise ValueError("df must be positive")
    if x <= 0:
        return 1.0
    return betainc(0.5 * df2, 0.5 * df1, df2 / (df2 + df1 * x))


def f_cdf(x, df1, df2):
    return 1.0 - f_sf(x, df1, df2)


# --- Standard normal ------------------------------------------------------------

def norm_cdf(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def norm_ppf(p):
    """Inverse of the standard normal CDF, by bisection on erfc (full double precision)."""
    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    if p == 0.5:
        return 0.0
    if p > 0.5:
        return -norm_ppf(1.0 - p)
    lo, hi = -40.0, 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid in (lo, hi):
            break
        if norm_cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
