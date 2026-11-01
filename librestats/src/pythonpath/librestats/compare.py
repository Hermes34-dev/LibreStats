"""Comparison tests: t-tests (one-sample, Welch, paired) and ANOVA (one-way,
Welch, repeated measures).

The t-tests follow Student (1908); the two-sample test is Welch's (1947) with
the Welch-Satterthwaite degrees of freedom, which does not assume equal
variances. Welch's ANOVA is the heteroscedastic F-test of Welch (1951).

The repeated-measures ANOVA is the univariate one-factor within-subjects
design (subjects as blocks). Sphericity is tested with Mauchly's (1940) W on
orthonormal contrasts; the Greenhouse-Geisser (1959) and Huynh-Feldt (1976)
epsilons correct the degrees of freedom. For a single group of subjects the
Lecoutre (1991) correction of the Huynh-Feldt epsilon changes nothing.
"""

import math
from dataclasses import dataclass, field

from .distributions import chi2_sf, f_sf, t_isf, t_sf
from .linalg import SingularMatrixError, cholesky

ALTERNATIVES = ("two-sided", "greater", "less")


class CompareError(ValueError):
    pass


@dataclass
class TTestResult:
    statistic: float
    df: float
    p: float
    estimate: float         # mean, or mean difference
    se: float
    ci: tuple               # (lower, upper); one-sided tests have an infinite bound
    effect: float           # Cohen's d (variant depends on the test)
    n: int
    notes: list = field(default_factory=list)


@dataclass
class AnovaResult:
    statistic: float
    df1: float
    df2: float
    p: float
    table: list = field(default_factory=list)       # (source, SS, df, MS, F or None, p or None)
    effects: list = field(default_factory=list)     # (label, value)
    sphericity: dict = field(default_factory=dict)  # repeated measures only
    corrections: list = field(default_factory=list)  # (name, epsilon, df1, df2, p)
    n: int = 0
    notes: list = field(default_factory=list)


def _clean(values):
    return [v for v in values if v is not None]


def _mean_var(x):
    n = len(x)
    mean = sum(x) / n
    return mean, sum((v - mean) ** 2 for v in x) / (n - 1)


def _t_test(t, df, estimate, se, alternative, alpha):
    """p-value and the (1 - alpha) confidence interval for a t statistic."""
    if alternative not in ALTERNATIVES:
        raise CompareError(f"Unknown alternative {alternative!r}.")
    if alternative == "two-sided":
        p = min(1.0, 2.0 * t_sf(abs(t), df))
        q = t_isf(alpha / 2, df)
        ci = (estimate - q * se, estimate + q * se)
    elif alternative == "greater":
        p = t_sf(t, df)
        ci = (estimate - t_isf(alpha, df) * se, math.inf)
    else:
        p = t_sf(-t, df)
        ci = (-math.inf, estimate + t_isf(alpha, df) * se)
    return p, ci


# --- t-tests ----------------------------------------------------------------------------

def one_sample_t(values, mu=0.0, alternative="two-sided", alpha=0.05):
    """One-sample t-test of H0: mean = mu. Effect size d = (mean - mu) / SD."""
    x = _clean(values)
    n = len(x)
    if n < 2:
        raise CompareError("The t-test needs at least 2 values.")
    mean, var = _mean_var(x)
    if var == 0:
        raise CompareError("All values are equal, so the standard error is zero.")
    sd = math.sqrt(var)
    se = sd / math.sqrt(n)
    t = (mean - mu) / se
    p, ci = _t_test(t, n - 1, mean, se, alternative, alpha)
    return TTestResult(t, n - 1, p, mean, se, ci, (mean - mu) / sd, n)


def welch_t(x1, x2, alternative="two-sided", alpha=0.05):
    """Welch's two-sample t-test of H0: mean1 = mean2 (difference = 1 - 2).
    Effect size: Cohen's d with the average of the two variances as standardiser."""
    x1, x2 = _clean(x1), _clean(x2)
    n1, n2 = len(x1), len(x2)
    if n1 < 2 or n2 < 2:
        raise CompareError("Each group needs at least 2 values.")
    m1, v1 = _mean_var(x1)
    m2, v2 = _mean_var(x2)
    a, b = v1 / n1, v2 / n2
    if a + b == 0:
        raise CompareError("Neither group varies, so the standard error is zero.")
    se = math.sqrt(a + b)
    df = (a + b) ** 2 / (a * a / (n1 - 1) + b * b / (n2 - 1))
    diff = m1 - m2
    t = diff / se
    p, ci = _t_test(t, df, diff, se, alternative, alpha)
    return TTestResult(t, df, p, diff, se, ci, diff / math.sqrt((v1 + v2) / 2), n1 + n2)


def paired_t(x1, x2, alternative="two-sided", alpha=0.05):
    """Paired t-test of H0: mean difference (1 - 2) = 0, on the pairs where both values
    are present. Effect size d_z = mean difference / SD of the differences."""
    if len(x1) != len(x2):
        raise CompareError("The two samples must have the same length.")
    d = [a - b for a, b in zip(x1, x2) if a is not None and b is not None]
    n = len(d)
    if n < 2:
        raise CompareError("The paired t-test needs at least 2 complete pairs.")
    mean, var = _mean_var(d)
    if var == 0:
        raise CompareError("All differences are equal, so the standard error is zero.")
    sd = math.sqrt(var)
    se = sd / math.sqrt(n)
    t = mean / se
    p, ci = _t_test(t, n - 1, mean, se, alternative, alpha)
    return TTestResult(t, n - 1, p, mean, se, ci, mean / sd, n)


# --- ANOVA ------------------------------------------------------------------------------

def _check_groups(groups, minimum=1):
    groups = [_clean(g) for g in groups]
    if len(groups) < 2:
        raise CompareError("The analysis needs at least 2 groups.")
    if any(len(g) < minimum for g in groups):
        raise CompareError(f"Every group needs at least {minimum} values.")
    return groups


def one_way_anova(groups):
    """Classical one-way ANOVA (Fisher). Effect sizes eta^2 and omega^2."""
    groups = _check_groups(groups)
    k = len(groups)
    n = sum(len(g) for g in groups)
    if n <= k:
        raise CompareError("The ANOVA needs more values than groups.")
    grand = sum(sum(g) for g in groups) / n
    means = [sum(g) / len(g) for g in groups]
    ss_b = sum(len(g) * (m - grand) ** 2 for g, m in zip(groups, means))
    ss_w = sum((v - m) ** 2 for g, m in zip(groups, means) for v in g)
    df1, df2 = k - 1, n - k
    if ss_w == 0:
        raise CompareError("No group varies, so the within-group variance is zero.")
    ms_b, ms_w = ss_b / df1, ss_w / df2
    f = ms_b / ms_w
    p = f_sf(f, df1, df2)
    ss_t = ss_b + ss_w
    table = [("Between groups", ss_b, df1, ms_b, f, p),
             ("Within groups", ss_w, df2, ms_w, None, None),
             ("Total", ss_t, n - 1, None, None, None)]
    effects = [("Eta squared (η²)", ss_b / ss_t),
               ("Omega squared (ω²)", (ss_b - df1 * ms_w) / (ss_t + ms_w))]
    return AnovaResult(f, df1, df2, p, table, effects, n=n)


def welch_anova(groups):
    """Welch's (1951) heteroscedastic one-way ANOVA."""
    groups = _check_groups(groups, 2)
    k = len(groups)
    stats = [(len(g),) + _mean_var(g) for g in groups]
    if any(v == 0 for _, _, v in stats):
        raise CompareError("Every group must vary (Welch's ANOVA weights groups by n / variance).")
    w = [n / v for n, _, v in stats]
    sw = sum(w)
    mw = sum(wi * m for wi, (_, m, _) in zip(w, stats)) / sw
    a = sum(wi * (m - mw) ** 2 for wi, (_, m, _) in zip(w, stats)) / (k - 1)
    h = sum((1 - wi / sw) ** 2 / (n - 1) for wi, (n, _, _) in zip(w, stats))
    f = a / (1 + 2 * (k - 2) * h / (k * k - 1))
    df1, df2 = k - 1, (k * k - 1) / (3 * h)
    return AnovaResult(f, df1, df2, f_sf(f, df1, df2), n=sum(n for n, _, _ in stats))


def _contrasts(k):
    """Orthonormal Helmert contrasts: k rows, k - 1 columns."""
    c = [[0.0] * (k - 1) for _ in range(k)]
    for j in range(1, k):
        s = math.sqrt(j * (j + 1))
        for i in range(j):
            c[i][j - 1] = 1 / s
        c[j][j - 1] = -j / s
    return c


def _sphericity(rows, n, k):
    """Mauchly's test and the Greenhouse-Geisser / Huynh-Feldt epsilons."""
    means = [sum(r[j] for r in rows) / n for j in range(k)]
    cov = [[sum((r[i] - means[i]) * (r[j] - means[j]) for r in rows) / (n - 1)
            for j in range(k)] for i in range(k)]
    c = _contrasts(k)
    p = k - 1
    cs = [[sum(c[a][i] * cov[a][b] for a in range(k)) for b in range(k)] for i in range(p)]
    sc = [[sum(cs[i][b] * c[b][j] for b in range(k)) for j in range(p)] for i in range(p)]
    tr = sum(sc[i][i] for i in range(p))
    tr2 = sum(sc[i][j] * sc[j][i] for i in range(p) for j in range(p))
    gg = tr * tr / (p * tr2)
    num, den = n * p * gg - 2, p * (n - 1 - p * gg)
    hf = 1.0 if den <= 0 else min(1.0, num / den)
    result = {"gg": gg, "hf": hf}
    if p >= 2 and n > p:
        try:
            log_det = 2 * sum(math.log(row[i]) for i, row in enumerate(cholesky(sc, 1e-14)))
        except SingularMatrixError:
            return result
        log_w = log_det - p * math.log(tr / p)
        chi2 = -(n - 1 - (2 * p * p + p + 2) / (6 * p)) * log_w
        df = p * (p + 1) // 2 - 1
        result.update(W=math.exp(log_w), chi2=chi2, df=df, p=chi2_sf(max(chi2, 0.0), df))
    return result


def rm_anova(rows):
    """One-factor repeated-measures ANOVA. rows: one list per subject, one value per
    condition; subjects with a missing value are dropped."""
    rows = [list(r) for r in rows if all(v is not None for v in r)]
    n = len(rows)
    k = len(rows[0]) if rows else 0
    if k < 2:
        raise CompareError("The repeated-measures ANOVA needs at least 2 conditions.")
    if n < 2:
        raise CompareError("The repeated-measures ANOVA needs at least 2 complete subjects.")
    grand = sum(sum(r) for r in rows) / (n * k)
    cond = [sum(r[j] for r in rows) / n for j in range(k)]
    subj = [sum(r) / k for r in rows]
    ss_t = sum((v - grand) ** 2 for r in rows for v in r)
    ss_c = n * sum((m - grand) ** 2 for m in cond)
    ss_s = k * sum((m - grand) ** 2 for m in subj)
    ss_e = ss_t - ss_c - ss_s
    df1, df2 = k - 1, (k - 1) * (n - 1)
    if ss_e <= 1e-12 * ss_t or ss_t == 0:
        raise CompareError("The error variance is zero: every subject changes by the same "
                           "amount between conditions.")
    ms_c, ms_e = ss_c / df1, ss_e / df2
    f = ms_c / ms_e
    p = f_sf(f, df1, df2)
    table = [("Conditions", ss_c, df1, ms_c, f, p),
             ("Subjects", ss_s, n - 1, ss_s / (n - 1), None, None),
             ("Error", ss_e, df2, ms_e, None, None),
             ("Total", ss_t, n * k - 1, None, None, None)]
    effects = [("Partial eta squared (η²p)", ss_c / (ss_c + ss_e)),
               ("Generalized eta squared (η²G)", ss_c / (ss_c + ss_s + ss_e))]
    sph = _sphericity(rows, n, k)
    corrections = [("Sphericity assumed", 1.0, df1, df2, p)]
    for name, eps in (("Greenhouse-Geisser", sph["gg"]), ("Huynh-Feldt", sph["hf"]),
                      ("Lower bound", 1 / (k - 1))):
        corrections.append((name, eps, eps * df1, eps * df2, f_sf(f, eps * df1, eps * df2)))
    notes = []
    if k == 2:
        notes.append("With 2 conditions sphericity always holds, and the F-test equals the "
                     "squared paired t-test.")
    elif "W" not in sph:
        notes.append("Mauchly's test needs more subjects than conditions minus one, and a "
                     "non-singular covariance matrix; it was not computed.")
    return AnovaResult(f, df1, df2, p, table, effects, sph, corrections, n, notes)
