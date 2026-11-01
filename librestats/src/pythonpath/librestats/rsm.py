"""Analysis of response surface designs: quadratic model, stationary point, optimum.

The model is fitted by least squares in coded units:
    y = b0 + sum bi xi + sum bii xi^2 + sum bij xi xj (+ block terms)
Blocks use sum-to-zero (effect) coding, so b0 and every prediction refer to an
average block.

- ANOVA by groups of terms (linear, square, two-factor interaction, blocks):
  adjusted sums of squares, i.e. the increase in the residual sum of squares
  when the group is removed from the full model. The residual is split into
  lack of fit and pure error (repeated points within a block).
- Predicted R^2 from PRESS = sum (e_i / (1 - h_ii))^2.
- Stationary point (Box & Wilson, 1951; Myers et al., 2016, ch. 6): with
  y = b0 + x'b + x'Bx, xs = -B^-1 b / 2, ys = b0 + xs'b / 2. Canonical
  analysis: the eigenvalues of B; all negative = maximum, all positive =
  minimum, mixed signs = saddle point; an eigenvalue near 0 = a ridge.
- Optimum within the experimental region (the box of the coded settings run):
  the maximum and minimum of a quadratic over a box lie at a stationary point
  of the function restricted to some face of the box, so every face (each
  factor at its lower bound, upper bound or free) is solved exactly.
"""

import math
from dataclasses import dataclass, field
from itertools import combinations, product

from .distributions import f_sf, t_sf
from .doe import LETTERS
from .linalg import SingularMatrixError, inverse_spd, matvec

MODELS = ("quadratic", "interaction", "linear")
RIDGE = 1e-6          # |eigenvalue| below this share of the largest: a ridge
MAX_FACES_K = 9       # 3^k faces


class RSMError(ValueError):
    pass


@dataclass
class Term:
    name: str                 # "A", "A^2", "AB", "Block 1"
    kind: str                 # linear | square | interaction | block
    factors: tuple            # factor indices
    coef: float = None
    se: float = None
    t: float = None
    p: float = None


@dataclass
class AnovaRow:
    source: str
    df: int
    ss: float
    f: float = None
    p: float = None


@dataclass
class Stationary:
    coded: list
    response: float
    eigenvalues: list
    eigenvectors: list        # columns, as lists
    nature: str               # maximum | minimum | saddle point
    ridge: bool
    inside: bool


@dataclass
class Fit:
    k: int
    model: str
    terms: list               # Term: constant first, then blocks, linear, square, interaction
    n: int
    missing: int
    df: int
    mse: float
    r2: float
    r2_adj: float
    r2_pred: float
    anova: list
    bounds: list              # (low, high) coded per factor: the experimental region
    stationary: Stationary = None
    maximum: tuple = None     # (coded point, predicted response)
    minimum: tuple = None
    notes: list = field(default_factory=list)

    def predict(self, x):
        """Prediction at coded point x for an average block."""
        return sum(t.coef * _value(t, x) for t in self.terms if t.kind != "block")


def _value(term, x):
    if term.kind == "constant":
        return 1.0
    if term.kind == "linear":
        return x[term.factors[0]]
    if term.kind == "square":
        return x[term.factors[0]] ** 2
    if term.kind == "interaction":
        return x[term.factors[0]] * x[term.factors[1]]
    return 0.0


def model_terms(k, model):
    if model not in MODELS:
        raise RSMError(f"Unknown model {model!r}.")
    terms = [Term("Constant", "constant", ())]
    terms += [Term(LETTERS[i], "linear", (i,)) for i in range(k)]
    if model == "quadratic":
        terms += [Term(LETTERS[i] + "^2", "square", (i,)) for i in range(k)]
    if model in ("quadratic", "interaction"):
        terms += [Term(LETTERS[i] + LETTERS[j], "interaction", (i, j))
                  for i, j in combinations(range(k), 2)]
    return terms


def _solve(a, b):
    """Gaussian elimination with partial pivoting; None if (nearly) singular."""
    n = len(b)
    m = [list(map(float, row)) + [float(v)] for row, v in zip(a, b)]
    scale = max((abs(v) for row in a for v in row), default=0.0) or 1.0
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        if abs(m[p][c]) < 1e-12 * scale:
            return None
        m[c], m[p] = m[p], m[c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            for j in range(c, n + 1):
                m[r][j] -= f * m[c][j]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (m[r][n] - sum(m[r][j] * x[j] for j in range(r + 1, n))) / m[r][r]
    return x


def eigen_symmetric(a, tol=1e-30, sweeps=100):
    """Eigenvalues and eigenvectors (columns) of a symmetric matrix (cyclic Jacobi)."""
    n = len(a)
    a = [list(map(float, row)) for row in a]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < tol * max(1.0, sum(a[i][i] ** 2 for i in range(n))):
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if a[p][q] == 0:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for r in range(n):
                    arp, arq = a[r][p], a[r][q]
                    a[r][p], a[r][q] = c * arp - s * arq, s * arp + c * arq
                for r in range(n):
                    apr, aqr = a[p][r], a[q][r]
                    a[p][r], a[q][r] = c * apr - s * aqr, s * apr + c * aqr
                for r in range(n):
                    vrp, vrq = v[r][p], v[r][q]
                    v[r][p], v[r][q] = c * vrp - s * vrq, s * vrp + c * vrq
    order = sorted(range(n), key=lambda i: a[i][i])
    return [a[i][i] for i in order], [[v[r][i] for r in range(n)] for i in order]


def quadratic_parts(fit):
    """(b0, b, B) of the fitted model in coded units: y = b0 + x'b + x'Bx."""
    k = fit.k
    b0 = next(t.coef for t in fit.terms if t.kind == "constant")
    b = [0.0] * k
    B = [[0.0] * k for _ in range(k)]
    for t in fit.terms:
        if t.kind == "linear":
            b[t.factors[0]] = t.coef
        elif t.kind == "square":
            B[t.factors[0]][t.factors[0]] = t.coef
        elif t.kind == "interaction":
            i, j = t.factors
            B[i][j] = B[j][i] = t.coef / 2
    return b0, b, B


def stationary_point(fit):
    b0, b, B = quadratic_parts(fit)
    values, vectors = eigen_symmetric(B)
    biggest = max(abs(v) for v in values)
    if biggest == 0:
        return None
    ridge = any(abs(v) < RIDGE * biggest for v in values)
    xs = _solve(B, [-v / 2 for v in b])
    if xs is None:
        return Stationary(None, None, values, vectors, "none (ridge)", True, False)
    ys = b0 + sum(x * v for x, v in zip(xs, b)) / 2
    if all(v < 0 for v in values):
        nature = "maximum"
    elif all(v > 0 for v in values):
        nature = "minimum"
    else:
        nature = "saddle point"
    inside = all(lo - 1e-9 <= x <= hi + 1e-9 for x, (lo, hi) in zip(xs, fit.bounds))
    return Stationary(xs, ys, values, vectors, nature, ridge, inside)


def box_optimum(fit, maximize=True):
    """Exact maximum (or minimum) of the fitted model over the coded box fit.bounds."""
    b0, b, B = quadratic_parts(fit)
    k = fit.k
    sign = 1 if maximize else -1
    best = None
    for state in product((0, 1, 2), repeat=k):           # 0 lower, 1 upper, 2 free
        x = [fit.bounds[i][state[i]] if state[i] < 2 else 0.0 for i in range(k)]
        free = [i for i in range(k) if state[i] == 2]
        if free:
            # gradient of the free coordinates: b_f + 2 B_ff x_f + 2 B_fF x_F = 0
            rhs = [-(b[i] + 2 * sum(B[i][j] * x[j] for j in range(k) if j not in free)) / 2
                   for i in free]
            sol = _solve([[B[i][j] for j in free] for i in free], rhs)
            if sol is None:
                continue
            for i, v in zip(free, sol):
                x[i] = v
            if any(not (fit.bounds[i][0] - 1e-12 <= x[i] <= fit.bounds[i][1] + 1e-12)
                   for i in free):
                continue
        y = fit.predict(x)
        if best is None or sign * y > sign * best[1]:
            best = (x, y)
    return best


def _least_squares(rows, y):
    p = len(rows[0])
    xtx = [[sum(r[a] * r[c] for r in rows) for c in range(p)] for a in range(p)]
    xty = [sum(r[a] * v for r, v in zip(rows, y)) for a in range(p)]
    inv = inverse_spd(xtx)
    beta = matvec(inv, xty)
    sse = sum((v - sum(a * c for a, c in zip(r, beta))) ** 2 for r, v in zip(rows, y))
    return beta, inv, sse


def fit_model(coded, blocks, y, model="quadratic"):
    """Fit the model to the runs with a response (y None = missing)."""
    k = len(coded[0])
    obs = [i for i, v in enumerate(y) if v is not None]
    missing = len(y) - len(obs)
    terms = model_terms(k, model)
    levels = sorted({blocks[i] for i in obs})
    block_terms = [Term(f"Block {b}", "block", (j,)) for j, b in enumerate(levels[:-1])]
    terms = terms[:1] + block_terms + terms[1:]

    def row(i):
        out = []
        for t in terms:
            if t.kind == "block":
                out.append(1.0 if blocks[i] == levels[t.factors[0]] else
                           -1.0 if blocks[i] == levels[-1] else 0.0)
            else:
                out.append(float(_value(t, coded[i])))
        return out

    rows = [row(i) for i in obs]
    yy = [y[i] for i in obs]
    n, p = len(rows), len(terms)
    if n <= p:
        raise RSMError(f"{n} runs with a response are too few for the {p} terms of the "
                       f"{model} model; at least {p + 1} are needed.")
    try:
        beta, inv, sse = _least_squares(rows, yy)
    except SingularMatrixError:
        hint = (" A two-level design with centre points cannot separate the square terms: "
                "choose the model with interactions (no squares), or add axial points."
                if model == "quadratic" else "")
        raise RSMError("The design cannot estimate all terms of the model." + hint) from None
    df = n - p
    mse = sse / df
    mean = sum(yy) / n
    sst = sum((v - mean) ** 2 for v in yy)
    for t, bcoef, j in zip(terms, beta, range(p)):
        t.coef = bcoef
        t.se = math.sqrt(mse * inv[j][j])
        t.t = bcoef / t.se if t.se > 0 else math.inf
        t.p = 2 * t_sf(abs(t.t), df)
    r2 = 1 - sse / sst if sst > 0 else None
    r2_adj = 1 - mse / (sst / (n - 1)) if sst > 0 else None
    press, leverage_one = 0.0, False
    for r, v in zip(rows, yy):
        h = sum(r[a] * sum(inv[a][c] * r[c] for c in range(p)) for a in range(p))
        if h > 1 - 1e-10:
            leverage_one = True
            break
        e = v - sum(a * c for a, c in zip(r, beta))
        press += (e / (1 - h)) ** 2
    r2_pred = None if leverage_one or sst <= 0 else 1 - press / sst

    # ANOVA by groups: adjusted SS = SSE without the group - SSE
    anova = []
    groups = [("Blocks", "block"), ("Linear", "linear"), ("Square", "square"),
              ("2-way interaction", "interaction")]
    model_df = p - 1
    model_ss = sst - sse
    anova.append(AnovaRow("Model", model_df, model_ss, (model_ss / model_df) / mse,
                          f_sf((model_ss / model_df) / mse, model_df, df)))
    for label, kind in groups:
        idx = [j for j, t in enumerate(terms) if t.kind == kind]
        if not idx:
            continue
        keep = [j for j in range(p) if j not in idx]
        _, _, sse_reduced = _least_squares([[r[j] for j in keep] for r in rows], yy)
        ss = sse_reduced - sse
        f = (ss / len(idx)) / mse
        anova.append(AnovaRow(label, len(idx), ss, f, f_sf(f, len(idx), df)))
    anova.append(AnovaRow("Residual error", df, sse))
    groups_pe = {}
    for i in obs:
        groups_pe.setdefault((blocks[i], tuple(coded[i])), []).append(y[i])
    pure_ss = sum(sum((v - sum(g) / len(g)) ** 2 for v in g) for g in groups_pe.values())
    pure_df = sum(len(g) - 1 for g in groups_pe.values())
    notes = []
    if pure_df and df > pure_df:
        lof_df, lof_ss = df - pure_df, sse - pure_ss
        f = (lof_ss / lof_df) / (pure_ss / pure_df) if pure_ss > 0 else None
        anova.append(AnovaRow("  Lack of fit", lof_df, lof_ss, f,
                              f_sf(f, lof_df, pure_df) if f is not None else None))
        anova.append(AnovaRow("  Pure error", pure_df, pure_ss))
    elif not pure_df:
        notes.append("No repeated points: lack of fit cannot be tested.")
    anova.append(AnovaRow("Total", n - 1, sst))
    if missing:
        notes.append(f"{missing} run(s) without a response are left out.")

    bounds = [(min(coded[i][j] for i in obs), max(coded[i][j] for i in obs))
              for j in range(k)]
    fit = Fit(k, model, terms, n, missing, df, mse, r2, r2_adj, r2_pred, anova, bounds,
              notes=notes)
    if model == "quadratic":
        fit.stationary = stationary_point(fit)
    if k <= MAX_FACES_K:
        fit.maximum = box_optimum(fit, True)
        fit.minimum = box_optimum(fit, False)
    return fit


def actual_equation(fit, centres, halves):
    """Coefficients in actual units, x = (X - c) / h: [(name, coefficient)] with
    names "Constant", "A", "A^2", "AB" (for an average block)."""
    b0, b, B = quadratic_parts(fit)
    k = fit.k
    c, h = centres, halves
    const = b0 + sum(-b[i] * c[i] / h[i] + B[i][i] * c[i] ** 2 / h[i] ** 2 for i in range(k))
    const += sum(2 * B[i][j] * c[i] * c[j] / (h[i] * h[j]) for i, j in combinations(range(k), 2))
    out = [("Constant", const)]
    for i in range(k):
        lin = b[i] / h[i] - 2 * B[i][i] * c[i] / h[i] ** 2
        lin -= sum(2 * B[i][j] * c[j] / (h[i] * h[j]) for j in range(k) if j != i)
        out.append((LETTERS[i], lin))
    kinds = {t.name for t in fit.terms}
    for i in range(k):
        if LETTERS[i] + "^2" in kinds:
            out.append((LETTERS[i] + "^2", B[i][i] / h[i] ** 2))
    for i, j in combinations(range(k), 2):
        if LETTERS[i] + LETTERS[j] in kinds:
            out.append((LETTERS[i] + LETTERS[j], 2 * B[i][j] / (h[i] * h[j])))
    return out


def _rounded_setting(value, half):
    """LibreStats' rounding of actual settings on a Design sheet (doe.Factor.setting)."""
    return float(f"{round(value, max(0, 3 - math.floor(math.log10(abs(half))))):.12g}")


def _exact_coding(xs, ys):
    """Centre and half-range from the rows at coded +-1 (Low/High of a circumscribed
    design, written exactly) or at the extreme coded levels (Low/High of an inscribed
    design), if every row's actual value then matches, exactly or as LibreStats rounds
    it. None if neither fits. When both fit (rounding can hide the difference), the
    pair with fewer decimals wins: Low and High are the numbers typed by the user,
    the other settings were rounded by LibreStats to more digits."""
    top = max(abs(x) for x in xs)
    found = []
    for level in (1.0, top):
        plus = [y for x, y in zip(xs, ys) if abs(x - level) < 1e-9]
        minus = [y for x, y in zip(xs, ys) if abs(x + level) < 1e-9]
        if not plus or not minus or len(set(plus)) > 1 or len(set(minus)) > 1:
            continue
        centre, half = (plus[0] + minus[0]) / 2, (plus[0] - minus[0]) / (2 * level)
        if half == 0:
            continue
        rounding_half = half * level        # the factor's (High - Low) / 2
        if all(abs(centre + half * x - y) <= 1e-9 * max(1.0, abs(y)) or
               _rounded_setting(centre + half * x, rounding_half) == y
               for x, y in zip(xs, ys)):
            found.append((_decimals(plus[0]) + _decimals(minus[0]), level != 1.0,
                          centre, half))
    return min(found)[2:] if found else None


def _decimals(v):
    text = f"{v:.12g}"
    return 0 if "." not in text or "e" in text else len(text.split(".")[1])


def coding(coded, actual):
    """Centre and half-range of each factor from the coded and actual columns
    (actual = centre + half * coded); also the largest misfit of the least-squares
    fallback. The actual settings on a Design sheet are rounded, so the exact Low/High
    rows are used when they can be identified."""
    k = len(coded[0])
    centres, halves, worst = [], [], 0.0
    for j in range(k):
        xs = [r[j] for r in coded]
        ys = [r[j] for r in actual]
        exact = _exact_coding(xs, ys)
        if exact:
            centres.append(exact[0])
            halves.append(exact[1])
            continue
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        sxx = sum((x - mx) ** 2 for x in xs)
        if sxx == 0:
            raise RSMError(f"Factor {LETTERS[j]} has only one level.")
        slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
        centre = my - slope * mx
        span = max(ys) - min(ys) or 1.0
        worst = max(worst, max(abs(centre + slope * x - y) for x, y in zip(xs, ys)) / span)
        centres.append(centre)
        halves.append(slope)
    return centres, halves, worst
