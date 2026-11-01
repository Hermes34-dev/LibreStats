"""Post hoc multiple comparisons after an ANOVA or a Kruskal-Wallis test.

Tukey's HSD (Tukey, 1953) compares all pairs of means with the pooled
within-group variance and the studentized range distribution; with unequal
group sizes it is the Tukey-Kramer procedure (Kramer, 1956). Games-Howell
(1976) does the same with each pair's own variances and Welch degrees of
freedom, so it does not assume equal variances. Dunnett's (1955) test compares
each group with one control group using the joint multivariate t distribution
of the comparisons. Dunn's (1964) test compares mean ranks after a
Kruskal-Wallis test, with the p-values adjusted for multiple testing.

All adjusted p-values and confidence intervals are simultaneous (family-wise)
for the set of comparisons made.
"""

import math
from dataclasses import dataclass, field

from .compare import ALTERNATIVES, CompareError, _check_groups, _mean_var
from .distributions import chi2_sf, dunnett_isf, dunnett_sf, tukey_isf, tukey_sf
from .transform import ranks

P_ADJUST = ("bonferroni", "sidak", "holm", "bh")


@dataclass
class Comparison:
    first: int              # group indices; estimate = first - second
    second: int
    estimate: float         # difference of means, or of mean ranks (Dunn)
    se: float
    statistic: float        # q (Tukey), t (Games-Howell, Dunnett), z (Dunn)
    df: float = None
    p: float = None         # adjusted p-value
    ci: tuple = None
    p_raw: float = None     # unadjusted p-value (Dunn)


@dataclass
class PosthocResult:
    comparisons: list
    params: list = field(default_factory=list)   # (label, value) for the report
    mean_ranks: list = None                      # Dunn only


def _pooled(groups):
    """Means, sizes, pooled within-group variance (MSE) and its df."""
    k = len(groups)
    n = sum(len(g) for g in groups)
    if n <= k:
        raise CompareError("The comparison needs more values than groups.")
    means = [sum(g) / len(g) for g in groups]
    ss_w = sum((v - m) ** 2 for g, m in zip(groups, means) for v in g)
    if ss_w == 0:
        raise CompareError("No group varies, so the within-group variance is zero.")
    return means, [len(g) for g in groups], ss_w / (n - k), n - k


def _pairs(k):
    return [(i, j) for i in range(k) for j in range(i + 1, k)]


def tukey_hsd(groups, alpha=0.05):
    """All pairwise comparisons with Tukey's HSD (Tukey-Kramer for unequal n)."""
    groups = _check_groups(groups)
    k = len(groups)
    means, ns, mse, df = _pooled(groups)
    q_crit = tukey_isf(alpha, k, df)
    out = []
    for i, j in _pairs(k):
        diff = means[i] - means[j]
        se = math.sqrt(mse * (1 / ns[i] + 1 / ns[j]))
        q = math.sqrt(2) * abs(diff) / se
        half = q_crit / math.sqrt(2) * se
        out.append(Comparison(i, j, diff, se, q, df, tukey_sf(q, k, df), (diff - half, diff + half)))
    return PosthocResult(out, [("Pooled variance (MSE)", mse), ("df", df),
                               ("Critical q", q_crit)])


def games_howell(groups, alpha=0.05):
    """All pairwise comparisons with the Games-Howell procedure (unequal variances)."""
    groups = _check_groups(groups, 2)
    k = len(groups)
    stats = [(len(g),) + _mean_var(g) for g in groups]
    if any(v == 0 for _, _, v in stats):
        raise CompareError("Every group must vary (Games-Howell uses each group's variance).")
    out = []
    for i, j in _pairs(k):
        (ni, mi, vi), (nj, mj, vj) = stats[i], stats[j]
        a, b = vi / ni, vj / nj
        se = math.sqrt(a + b)
        df = (a + b) ** 2 / (a * a / (ni - 1) + b * b / (nj - 1))
        diff = mi - mj
        t = diff / se
        half = tukey_isf(alpha, k, df) / math.sqrt(2) * se
        out.append(Comparison(i, j, diff, se, t, df, tukey_sf(math.sqrt(2) * abs(t), k, df),
                              (diff - half, diff + half)))
    return PosthocResult(out)


def dunnett(groups, control=0, alternative="two-sided", alpha=0.05):
    """Each group against the control group (index control), Dunnett's test."""
    if alternative not in ALTERNATIVES:
        raise CompareError(f"Unknown alternative {alternative!r}.")
    groups = _check_groups(groups)
    k = len(groups)
    if not 0 <= control < k:
        raise CompareError("The control group is not one of the groups.")
    means, ns, mse, df = _pooled(groups)
    others = [i for i in range(k) if i != control]
    n0 = ns[control]
    lambdas = [math.sqrt(ns[i] / (ns[i] + n0)) for i in others]
    two_sided = alternative == "two-sided"
    crit = dunnett_isf(alpha, df, lambdas, two_sided)
    out = []
    for i in others:
        diff = means[i] - means[control]
        se = math.sqrt(mse * (1 / ns[i] + 1 / n0))
        t = diff / se
        if two_sided:
            p = dunnett_sf(abs(t), df, lambdas)
            ci = (diff - crit * se, diff + crit * se)
        elif alternative == "greater":
            p = dunnett_sf(t, df, lambdas, False)
            ci = (diff - crit * se, math.inf)
        else:
            p = dunnett_sf(-t, df, lambdas, False)
            ci = (-math.inf, diff + crit * se)
        out.append(Comparison(i, control, diff, se, t, df, p, ci))
    return PosthocResult(out, [("Pooled variance (MSE)", mse), ("df", df),
                               ("Critical t", crit)])


# --- Rank-based: Kruskal-Wallis and Dunn ------------------------------------------------

def _ranked(groups):
    pooled = [v for g in groups for v in g]
    r = ranks(pooled).values
    out, start = [], 0
    for g in groups:
        out.append(r[start:start + len(g)])
        start += len(g)
    ties = {}
    for v in pooled:
        ties[v] = ties.get(v, 0) + 1
    return out, len(pooled), sum(t ** 3 - t for t in ties.values())


def kruskal_wallis(groups):
    """Kruskal-Wallis H (tie-corrected), df and p-value."""
    groups = _check_groups(groups)
    rk, n, tie_sum = _ranked(groups)
    if tie_sum == n ** 3 - n:
        raise CompareError("All values are equal.")
    h = 12 / (n * (n + 1)) * sum(sum(r) ** 2 / len(r) for r in rk) - 3 * (n + 1)
    h /= 1 - tie_sum / (n ** 3 - n)
    df = len(groups) - 1
    return h, df, chi2_sf(h, df)


def p_adjust(ps, method):
    """Adjust p-values for multiple testing: bonferroni, sidak, holm or bh
    (Benjamini-Hochberg false discovery rate)."""
    m = len(ps)
    if method == "bonferroni":
        return [min(1.0, m * p) for p in ps]
    if method == "sidak":
        return [-math.expm1(m * math.log1p(-p)) if p < 1 else 1.0 for p in ps]
    order = sorted(range(m), key=lambda i: ps[i])
    out = [0.0] * m
    if method == "holm":
        running = 0.0
        for rank, i in enumerate(order):
            running = max(running, min(1.0, (m - rank) * ps[i]))
            out[i] = running
        return out
    if method == "bh":
        running = 1.0
        for rank in reversed(range(m)):
            i = order[rank]
            running = min(running, m * ps[i] / (rank + 1))
            out[i] = running
        return out
    raise CompareError(f"Unknown p-value adjustment {method!r}.")


def dunn(groups, adjust="holm"):
    """Dunn's test of all pairs of mean ranks, with tie correction."""
    groups = _check_groups(groups)
    rk, n, tie_sum = _ranked(groups)
    var = n * (n + 1) / 12 - tie_sum / (12 * (n - 1))
    if var <= 0:
        raise CompareError("All values are equal.")
    mean_ranks = [sum(r) / len(r) for r in rk]
    out = []
    for i, j in _pairs(len(groups)):
        diff = mean_ranks[i] - mean_ranks[j]
        se = math.sqrt(var * (1 / len(rk[i]) + 1 / len(rk[j])))
        z = diff / se
        out.append(Comparison(i, j, diff, se, z, p_raw=math.erfc(abs(z) / math.sqrt(2))))
    for c, p in zip(out, p_adjust([c.p_raw for c in out], adjust)):
        c.p = p
    return PosthocResult(out, mean_ranks=mean_ranks)
