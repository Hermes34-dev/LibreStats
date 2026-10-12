"""Missing data: diagnostics (patterns, Little's MCAR test) and imputation (pure Python).

Data are a list of rows; each row is a list of floats with None for missing.

References
----------
Little, R. J. A. (1988). A test of missing completely at random for
    multivariate data with missing values. JASA, 83(404), 1198-1202.
Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from
    incomplete data via the EM algorithm. JRSS B, 39(1), 1-38.
Troyanskaya, O., et al. (2001). Missing value estimation methods for DNA
    microarrays. Bioinformatics, 17(6), 520-525.  (k-nearest-neighbour imputation)
van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
    imputation by chained equations in R. J. Stat. Softw., 45(3), 1-67.
Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
"""

import math
import random
import statistics
from dataclasses import dataclass, field

from . import distributions, linalg


class MissingDataError(ValueError):
    """The method cannot be applied to this data."""


def _observed(column):
    return [v for v in column if v is not None]


def _columns(rows, p):
    return [[r[j] for r in rows] for j in range(p)]


# --- Patterns -------------------------------------------------------------------

@dataclass
class Pattern:
    observed: tuple   # True where the variable is observed
    rows: list        # row indices with this pattern

    @property
    def count(self):
        return len(self.rows)


def patterns(rows):
    """Missing-data patterns, most frequent first."""
    groups = {}
    for i, r in enumerate(rows):
        groups.setdefault(tuple(v is not None for v in r), []).append(i)
    return sorted((Pattern(k, v) for k, v in groups.items()),
                  key=lambda pt: (-pt.count, [not o for o in pt.observed]))


# --- EM estimates of mean and covariance -------------------------------------------

@dataclass
class EMResult:
    mean: list
    cov: list          # maximum-likelihood covariance (divisor n)
    iterations: int
    n: int             # rows used (rows with at least one observed value)


def em_mvnorm(rows, tol=1e-10, max_iter=10000):
    """ML mean and covariance of multivariate normal data with missing values (EM)."""
    p = len(rows[0])
    rows = [r for r in rows if any(v is not None for v in r)]
    n = len(rows)
    cols = _columns(rows, p)
    obs_cols = [_observed(c) for c in cols]
    for j, c in enumerate(obs_cols):
        if len(c) < 2:
            raise MissingDataError(f"Variable {j + 1} has fewer than 2 observed values.")
    mean = [statistics.fmean(c) for c in obs_cols]
    cov = [[statistics.pvariance(obs_cols[i]) if i == j else 0.0 for j in range(p)]
           for i in range(p)]
    if any(cov[j][j] == 0 for j in range(p)):
        raise MissingDataError("A variable has zero variance.")
    groups = patterns(rows)

    for iteration in range(1, max_iter + 1):
        t1 = [0.0] * p
        t2 = [[0.0] * p for _ in range(p)]
        for pt in groups:
            o = [j for j in range(p) if pt.observed[j]]
            m = [j for j in range(p) if not pt.observed[j]]
            if m:
                try:
                    L = linalg.cholesky(linalg.submatrix(cov, o, o))
                except linalg.SingularMatrixError:
                    raise MissingDataError("The covariance matrix became singular "
                                           "(variables are collinear).") from None
                # B[k] = Soo^-1 Som[:, k]: regression of missing on observed
                B = [linalg.cho_solve(L, [cov[a][b] for a in o]) for b in m]
                cond = [[cov[a][b] - sum(cov[a][o[t]] * B[kb][t] for t in range(len(o)))
                         for kb, b in enumerate(m)] for a in m]
            for i in pt.rows:
                x = list(rows[i])
                if m:
                    dev = [x[a] - mean[a] for a in o]
                    for k, b in enumerate(m):
                        x[b] = mean[b] + sum(B[k][t] * dev[t] for t in range(len(o)))
                for a in range(p):
                    t1[a] += x[a]
                    for b in range(a, p):
                        t2[a][b] += x[a] * x[b]
            if m:
                for ka, a in enumerate(m):
                    for kb, b in enumerate(m):
                        if a <= b:
                            t2[a][b] += pt.count * cond[ka][kb]
        new_mean = [v / n for v in t1]
        new_cov = [[0.0] * p for _ in range(p)]
        for a in range(p):
            for b in range(a, p):
                new_cov[a][b] = new_cov[b][a] = t2[a][b] / n - new_mean[a] * new_mean[b]
        sd = [math.sqrt(max(new_cov[j][j], 1e-300)) for j in range(p)]
        change = max(
            max(abs(new_mean[j] - mean[j]) / sd[j] for j in range(p)),
            max(abs(new_cov[a][b] - cov[a][b]) / (sd[a] * sd[b])
                for a in range(p) for b in range(p)))
        mean, cov = new_mean, new_cov
        if change < tol:
            return EMResult(mean, cov, iteration, n)
    raise MissingDataError(f"EM did not converge in {max_iter} iterations.")


# --- Little's MCAR test ---------------------------------------------------------------

@dataclass
class LittleResult:
    statistic: float
    df: int
    p_value: float
    n_patterns: int
    n_rows: int
    em_iterations: int
    excluded_rows: int     # rows with every value missing


def little_mcar(rows):
    """Little's (1988) chi-square test of the null hypothesis that data are MCAR.

    Uses ML estimates of the mean and covariance from the EM algorithm, as in
    R's naniar::mcar_test. Assumes multivariate normality.
    """
    if not rows:
        raise MissingDataError("No data.")
    p = len(rows[0])
    if p < 2:
        raise MissingDataError("Little's test needs at least 2 variables.")
    used = [r for r in rows if any(v is not None for v in r)]
    if all(all(v is not None for v in r) for r in used):
        raise MissingDataError("There are no missing values, so Little's test does not apply.")
    em = em_mvnorm(used)
    d2, df = 0.0, -p
    groups = patterns(used)
    for pt in groups:
        o = [j for j in range(p) if pt.observed[j]]
        means = [statistics.fmean(used[i][j] for i in pt.rows) for j in o]
        diff = [means[t] - em.mean[j] for t, j in enumerate(o)]
        try:
            d2 += pt.count * linalg.quad_form_inv(linalg.submatrix(em.cov, o, o), diff)
        except linalg.SingularMatrixError:
            raise MissingDataError("The estimated covariance matrix is singular.") from None
        df += len(o)
    return LittleResult(d2, df, distributions.chi2_sf(d2, df), len(groups), len(used),
                        em.iterations, len(rows) - len(used))


# --- Single imputation ------------------------------------------------------------------

@dataclass
class Imputation:
    method: str
    data: list                                   # completed rows
    imputed: list = field(default_factory=list)  # (row, column) of each imputed cell
    summary: list = field(default_factory=list)  # [(variable index, label, value), ...]
    notes: list = field(default_factory=list)


def _missing_cells(rows):
    return [(i, j) for i, r in enumerate(rows) for j, v in enumerate(r) if v is None]


def impute_simple(rows, how="mean"):
    """Replace each missing value with its variable's mean or median."""
    p = len(rows[0])
    centre = statistics.fmean if how == "mean" else statistics.median
    fills = []
    for j, col in enumerate(_columns(rows, p)):
        obs = _observed(col)
        if not obs:
            raise MissingDataError(f"Variable {j + 1} has no observed values.")
        fills.append(centre(obs))
    data = [[fills[j] if v is None else v for j, v in enumerate(r)] for r in rows]
    res = Imputation(how.capitalize() + " imputation", data, _missing_cells(rows))
    res.summary = [(j, "Imputed value", fills[j]) for j in range(p)]
    res.notes.append(f"{how.capitalize()} imputation shrinks variances and correlations; "
                     "use it only when few values are missing.")
    return res


def impute_knn(rows, k=5, standardize=True):
    """k-nearest-neighbour imputation.

    For each missing value, the k nearest rows that have that variable observed
    (donors) are found using the Euclidean distance over the variables observed
    in both rows, scaled up for the variables not compared (as scikit-learn's
    nan_euclidean distance). The value is the mean of the donors' values.
    With `standardize`, distances use z-scores so variables count equally.
    """
    p = len(rows[0])
    cols = _columns(rows, p)
    centre, scale = [], []
    for j, col in enumerate(cols):
        obs = _observed(col)
        if not obs:
            raise MissingDataError(f"Variable {j + 1} has no observed values.")
        if standardize:
            sd = statistics.stdev(obs) if len(obs) > 1 else 0.0
            centre.append(statistics.fmean(obs))
            scale.append(sd if sd > 0 else 1.0)
        else:
            centre.append(0.0)
            scale.append(1.0)
    z = [[None if v is None else (v - centre[j]) / scale[j] for j, v in enumerate(r)]
         for r in rows]

    def distance(a, b):
        shared = [(x - y) ** 2 for x, y in zip(a, b) if x is not None and y is not None]
        if not shared:
            return None
        return math.sqrt(sum(shared) * p / len(shared))

    data = [list(r) for r in rows]
    fallback = 0
    for i, j in _missing_cells(rows):
        dists = []
        for d_row, other in enumerate(z):
            if other[j] is None:
                continue
            dist = distance(z[i], other)
            if dist is not None:
                dists.append((dist, d_row))
        dists.sort()
        donors = dists[:k]
        if donors:
            data[i][j] = statistics.fmean(rows[d][j] for _, d in donors)
        else:
            data[i][j] = statistics.fmean(_observed(cols[j]))
            fallback += 1
    res = Imputation(f"k-nearest neighbours (k = {k})", data, _missing_cells(rows))
    res.notes.append("Single imputation: standard errors computed from the completed data "
                     "will be too small.")
    if fallback:
        res.notes.append(f"{fallback} value(s) had no comparable donor row and were "
                         "replaced by the variable mean.")
    return res


# --- Multiple imputation by chained equations (MICE) ------------------------------------

def _bayes_regression_draw(X, y, rng):
    """Draw (beta_hat, beta_star, sigma_star) as mice.impute.norm.draw (Rubin 1987)."""
    q = len(X[0])
    n = len(y)
    xtx = [[sum(r[a] * r[b] for r in X) for b in range(q)] for a in range(q)]
    for a in range(q):  # small ridge penalty, as in mice, for numerical stability
        xtx[a][a] += 1e-5 * xtx[a][a] if xtx[a][a] > 0 else 1e-5
    xty = [sum(r[a] * yi for r, yi in zip(X, y)) for a in range(q)]
    L = linalg.cholesky(xtx)
    beta = linalg.cho_solve(L, xty)
    sse = sum((yi - sum(b * x for b, x in zip(beta, r))) ** 2 for r, yi in zip(X, y))
    g = rng.gammavariate((n - q) / 2.0, 2.0)  # chi-square with n - q df
    sigma = math.sqrt(sse / g)
    # beta* = beta + sigma * chol(V) z with V = (X'X)^-1. If X'X = L L^T then
    # V = L^-T L^-1, so L^-T z has covariance V.
    z = [rng.gauss(0.0, 1.0) for _ in range(q)]
    shift = linalg.backward_solve(L, z)
    beta_star = [b + sigma * s for b, s in zip(beta, shift)]
    return beta, beta_star, sigma


@dataclass
class MiceResult:
    datasets: list           # m completed copies of the data
    imputed: list            # (row, column) of each imputed cell
    method: str
    m: int
    iterations: int
    seed: object
    notes: list = field(default_factory=list)


def mice(rows, m=5, iterations=10, method="pmm", donors=5, seed=None):
    """Multiple imputation by chained equations for numeric variables.

    Each incomplete variable is imputed in turn from all other variables by
    Bayesian linear regression; method "norm" adds a residual draw, method
    "pmm" (predictive mean matching, mice's default) copies the observed
    value of one of the `donors` closest predicted cases. m independent chains
    of `iterations` iterations give m completed datasets.
    """
    if method not in ("pmm", "norm"):
        raise MissingDataError(f"Unknown MICE method {method!r}.")
    p = len(rows[0])
    if p < 2:
        raise MissingDataError("MICE needs at least 2 variables.")
    rows = [list(r) for r in rows]
    cols = _columns(rows, p)
    targets = [j for j in range(p) if any(v is None for v in cols[j])]
    if not targets:
        raise MissingDataError("There are no missing values to impute.")
    for j in targets:
        n_obs = len(_observed(cols[j]))
        if n_obs < p + 2:
            raise MissingDataError(f"Variable {j + 1} has only {n_obs} observed values; "
                                   f"MICE needs at least {p + 2} here.")
    missing = _missing_cells(rows)
    rng = random.Random(seed)
    datasets = []
    for _ in range(m):
        data = [list(r) for r in rows]
        for i, j in missing:  # start from random draws of observed values
            data[i][j] = rng.choice(_observed(cols[j]))
        for _ in range(iterations):
            for j in targets:
                obs_rows = [i for i in range(len(rows)) if rows[i][j] is not None]
                mis_rows = [i for i in range(len(rows)) if rows[i][j] is None]

                def design(i):
                    return [1.0] + [data[i][c] for c in range(p) if c != j]

                X_obs = [design(i) for i in obs_rows]
                y_obs = [rows[i][j] for i in obs_rows]
                try:
                    beta, beta_star, sigma = _bayes_regression_draw(X_obs, y_obs, rng)
                except linalg.SingularMatrixError:
                    raise MissingDataError("Regression failed: predictors are collinear "
                                           "or a variable is constant.") from None
                X_mis = [design(i) for i in mis_rows]
                pred_mis = [sum(b * x for b, x in zip(beta_star, r)) for r in X_mis]
                if method == "norm":
                    for i, pred in zip(mis_rows, pred_mis):
                        data[i][j] = pred + rng.gauss(0.0, sigma)
                else:
                    pred_obs = [sum(b * x for b, x in zip(beta, r)) for r in X_obs]
                    for i, pred in zip(mis_rows, pred_mis):
                        nearest = sorted(range(len(obs_rows)),
                                         key=lambda t: (abs(pred_obs[t] - pred), rng.random()))
                        data[i][j] = y_obs[rng.choice(nearest[:donors])]
        datasets.append(data)
    res = MiceResult(datasets, missing, method, m, iterations, seed)
    res.notes.append("Analyse each completed dataset separately and pool the results with "
                     "Rubin's rules; do not average the imputations into one dataset.")
    return res
