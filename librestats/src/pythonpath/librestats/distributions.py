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


# --- Studentized range and Dunnett's many-to-one t ----------------------------------------
#
# Both are mixtures over s = sqrt(chi2_df / df) of a normal-theory probability:
# P = integral f_df(s) U(s) ds. The outer integral is taken over t = log s, where
# the integrand is smooth and close to Gaussian, by composite Gauss-Legendre
# quadrature, and normalised by the integral of f_df on the same nodes. The
# inner integrals give the upper tail directly (not as 1 - CDF), so that small
# p-values keep their relative accuracy.

_GL_X = (-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472,
         -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244,
         0.8650633666889845, 0.9739065285171717)
_GL_W = (0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963,
         0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2190863625159820,
         0.1494513491505806, 0.0666713443086881)
_SQRT2 = math.sqrt(2.0)
_INV_SQRT_2PI = 1.0 / math.sqrt(2.0 * math.pi)


def _gl_nodes(a, b, width):
    """Composite 10-point Gauss-Legendre nodes and weights on [a, b]."""
    m = max(1, math.ceil((b - a) / width))
    h = (b - a) / m
    for j in range(m):
        mid = a + (j + 0.5) * h
        for x, w in zip(_GL_X, _GL_W):
            yield mid + 0.5 * h * x, 0.5 * h * w


def _norm_sf(x):
    return 0.5 * math.erfc(x / _SQRT2)


def _scale_mixture(tail, df, decay):
    """integral of f_df(s) tail(s) ds, s = sqrt(chi2_df / df); tail(s) ~ exp(-decay s^2)."""
    if df > 1e9:
        return tail(1.0)
    sigma = 1.0 / math.sqrt(2.0 * df)
    t_peak = 0.5 * math.log(df / (df + 2.0 * decay))
    lo = t_peak - max(9.0 * sigma, 45.0 / df)
    hi = max(9.0 * sigma, 0.5 * math.log(1.0 + 100.0 / df))
    total = norm = 0.0
    for t, w in _gl_nodes(lo, hi, min(1.0, 2.0 * sigma)):
        g = w * math.exp(df * t - 0.5 * df * (math.exp(2.0 * t) - 1.0))  # f(s) ds, unnormalised
        if g > 0.0:
            norm += g
            total += g * tail(math.exp(t))
    return total / norm


def _range_sf_normal(w, k):
    """P(range of k standard normals > w)."""
    if w <= 0.0:
        return 1.0
    m = k - 1
    total = 0.0
    centre = -0.5 * w
    for z, wt in _gl_nodes(centre - 9.0, min(9.0, centre + 9.0), 1.0):
        a = _norm_sf(z)
        c = _norm_sf(z + w)
        if c <= 0.0:
            continue
        # A^m - (A - c)^m, without cancellation
        inner = 1.0 if c >= a else -math.expm1(m * math.log1p(-c / a))
        total += wt * math.exp(-0.5 * z * z) * a ** m * inner
    return min(1.0, k * _INV_SQRT_2PI * total)


def tukey_sf(q, k, df):
    """Upper tail P(Q > q) of the studentized range for k means and df degrees of freedom."""
    if k < 2 or df <= 0:
        raise ValueError("k must be at least 2 and df positive")
    if q <= 0.0:
        return 1.0
    return min(1.0, _scale_mixture(lambda s: _range_sf_normal(q * s, k), df, q * q / 4.0))


def _invert_sf(sf, p, lo, hi):
    """x in [lo, hi] with sf(x) = p for a decreasing sf: Illinois method on log sf."""
    target = math.log(p)
    flo, fhi = math.log(sf(lo)) - target, math.log(sf(hi)) - target
    while fhi > 0:
        lo, flo, hi = hi, fhi, 2.0 * hi
        fhi = math.log(sf(hi)) - target
    side = 0
    for _ in range(100):
        x = (lo * fhi - hi * flo) / (fhi - flo)
        fx = math.log(sf(x)) - target
        if abs(fx) < 1e-12 or hi - lo < 1e-12 * hi:
            return x
        if fx > 0:
            lo, flo = x, fx
            if side == 1:
                fhi /= 2.0
            side = 1
        else:
            hi, fhi = x, fx
            if side == -1:
                flo /= 2.0
            side = -1
    return x


def tukey_isf(p, k, df):
    """q with P(Q > q) = p for the studentized range."""
    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    pairs = k * (k - 1) / 2
    lo = _SQRT2 * t_isf(min(0.5, p / 2.0), df)          # one comparison
    hi = _SQRT2 * t_isf(p / (2.0 * pairs), df) * 1.01  # Bonferroni bound
    return _invert_sf(lambda q: tukey_sf(q, k, df), p, max(lo, 1e-6), hi)


def _dunnett_tail_normal(c, lambdas, two_sided):
    """P(max_j |X_j| > c) (or max_j X_j > c), X_j standard normal with corr lambda_i lambda_j."""
    roots = [math.sqrt(1.0 - lam * lam) for lam in lambdas]
    reach = 9.0 + c * max(lambdas)
    total = 0.0
    for z, wt in _gl_nodes(-reach, reach, 1.0):
        log_all_inside = 0.0
        for lam, r in zip(lambdas, roots):
            e = _norm_sf((c - lam * z) / r)
            if two_sided:
                e += _norm_sf((c + lam * z) / r)
            log_all_inside += math.log1p(-min(e, 1.0)) if e < 1.0 else -math.inf
        total += wt * math.exp(-0.5 * z * z) * -math.expm1(log_all_inside)
    return min(1.0, _INV_SQRT_2PI * total)


def dunnett_sf(t, df, lambdas, two_sided=True):
    """P(max_j |T_j| > t) (two-sided) or P(max_j T_j > t) for Dunnett's many-to-one
    t statistics, with correlations lambda_i * lambda_j, lambda_i = sqrt(n_i / (n_i + n_0))."""
    if df <= 0:
        raise ValueError("df must be positive")
    if t <= 0.0 and two_sided:
        return 1.0
    if t <= -40.0:
        return 1.0
    return min(1.0, _scale_mixture(lambda s: _dunnett_tail_normal(t * s, lambdas, two_sided),
                                   df, max(t, 0.0) ** 2 / 2.0))


def dunnett_isf(p, df, lambdas, two_sided=True):
    """Critical value t with dunnett_sf(t) = p."""
    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    m = len(lambdas)
    side = 2.0 if two_sided else 1.0
    lo = t_isf(min(0.5, p / side), df)
    hi = t_isf(p / (side * m), df) * 1.01  # Bonferroni bound
    return _invert_sf(lambda t: dunnett_sf(t, df, lambdas, two_sided), p, max(lo, 1e-6), hi)
