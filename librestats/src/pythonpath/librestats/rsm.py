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
- R^2 ceiling = 1 - SS(pure error) / SS(total): the largest R^2 any model can reach,
  since no model explains the scatter of repeated runs (Draper & Smith, 1998).
- Per-run diagnostics with Cook's (1977) distance: see influence.py.
- Stationary point (Box & Wilson, 1951; Myers et al., 2016, ch. 6): with
  y = b0 + x'b + x'Bx, xs = -B^-1 b / 2, ys = b0 + xs'b / 2. Canonical
  analysis: the eigenvalues of B; all negative = maximum, all positive =
  minimum, mixed signs = saddle point; an eigenvalue near 0 = a ridge.
- Optimum within the experimental region (the box of the coded settings run):
  the maximum and minimum of a quadratic over a box lie at a stationary point
  of the function restricted to some face of the box, so every face (each
  factor at its lower bound, upper bound or free) is solved exactly.
- Reduced cubic: the full cubic (quadratic + A^3, A^2B, ABC terms) with the
  aliased terms removed. The columns of the model matrix are orthogonalised in
  hierarchical order (constant, blocks, linear, squares, two-factor
  interactions, then cubic terms: ABC-type, A^2B-type, A^3-type), i.e. a QR
  decomposition in a fixed column order: a column that is (numerically) a linear
  combination of the columns before it is aliased with them, is left out, and
  that combination is reported. A fixed order rather than column pivoting
  keeps the model hierarchical (a cubic term never displaces a lower-order
  one). Which member of an aliased set is kept does not change the fitted
  values, residuals, R^2 or ANOVA, only how the coefficients are labelled. The
  optimum of a cubic over the box is found numerically: a grid search, then
  coordinate-wise refinement (along one coordinate the model is at most cubic,
  so each step is solved exactly).
"""

import math
from dataclasses import dataclass, field
from itertools import combinations, combinations_with_replacement, product

from .distributions import f_sf, t_sf
from .doe import LETTERS
from .influence import diagnostics
from .linalg import SingularMatrixError, inverse_spd, matvec

MODELS = ("quadratic", "interaction", "linear", "cubic")
ALIAS_TOL = 1e-8      # residual norm / column norm below this: the column is aliased
GRID_POINTS = 20000   # grid search of the cubic optimum: about this many points
MAX_MAPS = 3          # response maps / 3D surfaces drawn at most
RIDGE = 1e-6          # |eigenvalue| below this share of the largest: a ridge
MAX_FACES_K = 9       # 3^k faces


class RSMError(ValueError):
    pass


@dataclass
class Term:
    name: str                 # "A", "A^2", "AB", "Block 1"
    kind: str                 # linear | square | interaction | cubic | block
    factors: tuple            # factor indices (cubic: with repeats, e.g. (0, 0, 1) = A^2B)
    coef: float = None
    se: float = None
    t: float = None
    p: float = None
    vif: float = None         # variance inflation factor (not for the constant)


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
    r2_max: float = None      # R^2 ceiling from pure error; None without repeated runs
    observations: list = field(default_factory=list)   # run indices that were fitted
    diagnostics: list = field(default_factory=list)    # influence.Diagnostic per fitted run
    aliases: list = field(default_factory=list)   # reduced cubic: Alias of each removed term
    start_terms: list = None  # names of all terms of the (alias-reduced) model, before keep
    lof_p: float = None       # lack-of-fit p (None if it cannot be tested)
    aicc: float = None
    full_terms: int = None    # reduced cubic: number of terms of the full cubic (with blocks)

    def predict(self, x):
        """Prediction at coded point x for an average block."""
        return sum(t.coef * _value(t, x) for t in self.terms if t.kind != "block")


@dataclass
class Alias:
    """A term left out of the model: term = sum of coef x kept term (coded units)."""
    term: str
    combination: list         # [(coefficient, term name)]

    def text(self):
        names = ["the constant" if n == "Constant" else n for _, n in self.combination]
        if not names:
            return f"{self.term} is zero in every run (it cannot be estimated)"
        if len(names) == 1 and abs(abs(self.combination[0][0]) - 1) < 1e-9:
            sign = "" if self.combination[0][0] > 0 else "minus "
            return f"{self.term} is aliased with {sign}{names[0]}"
        parts = []
        for i, (c, n) in enumerate(self.combination):
            mag = "" if abs(abs(c) - 1) < 1e-9 else f"{abs(c):.4g} "
            parts.append(("-" if c < 0 else "") * (i == 0) +
                         ("" if i == 0 else " - " if c < 0 else " + ") + mag + n)
        listed = ", ".join(names[:-1]) + " and " + names[-1]
        return f"{self.term} is aliased with {listed}: {self.term} = {''.join(parts)}"


def alias_marks(aliases):
    """{kept term: note} for the terms that stand in for an aliased term of the same
    order: 'aliased with AC^2' when the removed term equals it (up to sign), otherwise
    'part of the A^3 alias' (the removed term is a combination of several kept terms).
    Lower-order terms in a combination (A in A^3 = 2 A - AB^2) are not marked."""
    marks = {}
    for a in aliases:
        order = len(name_monomial(a.term))
        same = [(c, n) for c, n in a.combination
                if n != "Constant" and len(name_monomial(n)) == order]
        for c, n in same:
            if len(a.combination) == 1:
                note = f"aliased with {'-' if c < 0 else ''}{a.term}"
            else:
                note = f"part of the {a.term} alias"
            marks.setdefault(n, []).append(note)
    return {n: "; ".join(v) for n, v in marks.items()}


def corners_not_run(bounds, points, tol=1e-9):
    """Corners of the box (coded bounds) that no run is at."""
    run = [tuple(p) for p in points]
    out = []
    for corner in product(*bounds):
        if not any(all(abs(a - b) <= tol for a, b in zip(corner, p)) for p in run):
            out.append(corner)
    return out


def monomial(term):
    """The term as a sorted tuple of factor indices with repeats, e.g. A^2B = (0, 0, 1)."""
    if term.kind == "square":
        return (term.factors[0],) * 2
    return tuple(sorted(term.factors))


def monomial_name(m):
    """(0, 0, 1) -> 'A^2B'; () -> 'Constant'."""
    if not m:
        return "Constant"
    out = ""
    for i in sorted(set(m)):
        n = m.count(i)
        out += LETTERS[i] + (f"^{n}" if n > 1 else "")
    return out


def name_monomial(name):
    """'A^2B' -> (0, 0, 1); 'Constant' -> ()."""
    if name == "Constant":
        return ()
    out, i = [], 0
    while i < len(name):
        f = LETTERS.index(name[i])
        n = 1
        if name[i + 1:i + 2] == "^":
            n, i = int(name[i + 2]), i + 2
        out += [f] * n
        i += 1
    return tuple(sorted(out))


def _value(term, x):
    if term.kind == "cubic":
        return math.prod(x[i] for i in term.factors)
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
    if model == "cubic":
        terms += [Term(LETTERS[i] + "^2", "square", (i,)) for i in range(k)]
    if model in ("quadratic", "interaction", "cubic"):
        terms += [Term(LETTERS[i] + LETTERS[j], "interaction", (i, j))
                  for i, j in combinations(range(k), 2)]
    if model == "cubic":
        cubes = [m for m in combinations_with_replacement(range(k), 3)]
        cubes.sort(key=lambda m: (len(set(m)) != 3, len(set(m)) == 1, m))  # ABC, A^2B, A^3
        terms += [Term(monomial_name(m), "cubic", m) for m in cubes]
    return terms


def independent_columns(rows, names, tol=ALIAS_TOL):
    """Hierarchical rank check of a model matrix: Gram-Schmidt (twice, for accuracy) on
    the columns in the given order. Returns the indices of the columns kept and an
    Alias for each column that is a linear combination of the columns before it."""
    n, p = len(rows), len(rows[0])
    cols = [[r[j] for r in rows] for j in range(p)]
    q, rmat, kept, aliases = [], [], [], []      # rmat[j]: column j of R (kept columns)
    for j, v in enumerate(cols):
        norm = math.sqrt(sum(a * a for a in v))
        w, coefs = list(v), [0.0] * len(q)
        for _ in range(2):
            for t, qt in enumerate(q):
                c = sum(a * b for a, b in zip(qt, w))
                coefs[t] += c
                w = [a - c * b for a, b in zip(w, qt)]
        rest = math.sqrt(sum(a * a for a in w))
        if norm == 0 or rest <= tol * norm:
            # v = Q coefs = K a with K the kept columns, K = Q R: solve R a = coefs
            m = len(q)
            a = [0.0] * m
            for i in range(m - 1, -1, -1):
                a[i] = (coefs[i] - sum(rmat[c][i] * a[c] for c in range(i + 1, m))) / \
                    rmat[i][i]
            scale = max((abs(x) for x in a), default=0.0)
            combo = [(x, names[kept[i]]) for i, x in enumerate(a)
                     if abs(x) > 1e-9 * max(1.0, scale)]
            aliases.append(Alias(names[j], [(_tidy(x), nm) for x, nm in combo]))
            continue
        q.append([a / rest for a in w])
        rmat.append(coefs + [rest])
        kept.append(j)
    return kept, aliases


def _tidy(x):
    """Round a coefficient that is a whole number up to rounding error."""
    r = round(x)
    return float(r) if abs(x - r) < 1e-9 * max(1.0, abs(x)) else x


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


def fit_model(coded, blocks, y, model="quadratic", keep=None):
    """Fit the model to the runs with a response (y None = missing). keep: names of the
    terms to keep from the (alias-reduced) model, e.g. after model reduction; the
    constant and the block terms are always kept. None: all terms."""
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
    aliases, full_terms = [], None
    if model == "cubic":
        full_terms = len(terms)
        kept, aliases = independent_columns(rows, [t.name for t in terms])
        terms = [terms[j] for j in kept]
        rows = [[r[j] for j in kept] for r in rows]
    start_terms = [t.name for t in terms]
    if keep is not None:
        keep = set(keep)
        unknown = keep - set(start_terms)
        if unknown:
            raise RSMError("Not in the model: " + ", ".join(sorted(unknown)) + ".")
        idx = [j for j, t in enumerate(terms)
               if t.kind in ("constant", "block") or t.name in keep]
        terms = [terms[j] for j in idx]
        rows = [[r[j] for j in idx] for r in rows]
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
        if t.kind != "constant":
            col = [r[j] for r in rows]
            m = sum(col) / n
            t.vif = inv[j][j] * sum((v - m) ** 2 for v in col)
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
              ("2-way interaction", "interaction"), ("Cubic", "cubic")]
    model_df = p - 1
    model_ss = sst - sse
    if model_df:
        anova.append(AnovaRow("Model", model_df, model_ss, (model_ss / model_df) / mse,
                              f_sf((model_ss / model_df) / mse, model_df, df)))
    else:                                     # constant only (e.g. after model reduction)
        anova.append(AnovaRow("Model", 0, 0.0))
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
    lof = next((r for r in anova if r.source.strip() == "Lack of fit"), None)
    fit = Fit(k, model, terms, n, missing, df, mse, r2, r2_adj, r2_pred, anova, bounds,
              notes=notes, observations=obs,
              r2_max=1 - pure_ss / sst if pure_df and sst > 0 else None,
              aliases=aliases, full_terms=full_terms,
              diagnostics=diagnostics(rows, yy, beta, inv, sse),
              start_terms=start_terms, lof_p=lof.p if lof else None,
              aicc=aicc(sse, n, p))
    cubic = any(t.kind == "cubic" for t in terms)
    if not cubic and any(t.kind == "square" for t in terms):
        fit.stationary = stationary_point(fit)
    if cubic:
        fit.maximum = numeric_optimum(fit, True)
        fit.minimum = numeric_optimum(fit, False)
    elif k <= MAX_FACES_K:
        fit.maximum = box_optimum(fit, True)
        fit.minimum = box_optimum(fit, False)
    return fit


def aicc(sse, n, p):
    """AICc = n ln(SSE/n) + 2p + 2p(p+1)/(n-p-1), p = model terms including the constant
    (and block terms); None when n - p - 1 <= 0 or SSE = 0."""
    if n - p - 1 <= 0 or sse <= 0:
        return None
    return n * math.log(sse / n) + 2 * p + 2 * p * (p + 1) / (n - p - 1)


def map_pairs(fit, max_pairs=MAX_MAPS, points=21):
    """Factor pairs (i, j) for response maps: all pairs if there are at most max_pairs,
    otherwise the max_pairs with the largest range of predicted values over their map
    (the other factors at their centre), in factor order. Returns (pairs, all pairs)."""
    every = list(combinations(range(fit.k), 2))
    if len(every) <= max_pairs:
        return every, every

    def spread(pair):
        i, j = pair
        (x0, x1), (y0, y1) = fit.bounds[i], fit.bounds[j]
        values = []
        for a in range(points):
            for b in range(points):
                x = [0.0] * fit.k
                x[i] = x0 + (x1 - x0) * a / (points - 1)
                x[j] = y0 + (y1 - y0) * b / (points - 1)
                values.append(fit.predict(x))
        return max(values) - min(values)

    ranked = sorted(every, key=lambda pr: -spread(pr))
    return sorted(ranked[:max_pairs]), every


def _cubic_on_interval(f, lo, hi):
    """Points where f (a polynomial of degree <= 3 on [lo, hi]) can have its extremes:
    the ends and the roots of f' (f is interpolated from 4 points)."""
    if hi - lo <= 0:
        return [lo]
    ts = [lo + (hi - lo) * i / 3 for i in range(4)]
    vs = [f(t) for t in ts]
    # Newton divided differences -> power basis in u = (t - lo) / (hi - lo) * 3
    d1 = [vs[i + 1] - vs[i] for i in range(3)]
    d2 = [(d1[i + 1] - d1[i]) / 2 for i in range(2)]
    d3 = (d2[1] - d2[0]) / 3
    # f = v0 + d1[0] u + d2[0] u (u - 1) + d3 u (u - 1) (u - 2)
    a3 = d3
    a2 = d2[0] - 3 * d3
    a1 = d1[0] - d2[0] + 2 * d3
    # f'(u) = 3 a3 u^2 + 2 a2 u + a1
    roots = []
    if abs(a3) > 1e-14 * (abs(a2) + abs(a1) + 1e-300):
        disc = (2 * a2) ** 2 - 12 * a3 * a1
        if disc >= 0:
            sq = math.sqrt(disc)
            roots = [(-2 * a2 + sq) / (6 * a3), (-2 * a2 - sq) / (6 * a3)]
    elif abs(a2) > 0:
        roots = [-a1 / (2 * a2)]
    out = [lo, hi]
    for u in roots:
        t = lo + u * (hi - lo) / 3
        if lo < t < hi:
            out.append(t)
    return out


def numeric_optimum(fit, maximize=True, points=GRID_POINTS):
    """Maximum (or minimum) of any fitted model of degree <= 3 per factor over the
    coded box: the best points of a grid, each improved by coordinate-wise exact
    one-dimensional optimisation until no coordinate changes."""
    k, sign = fit.k, 1 if maximize else -1
    m = max(3, int(points ** (1 / k)))
    axes = [[lo + (hi - lo) * i / (m - 1) for i in range(m)] for lo, hi in fit.bounds]
    scored = sorted(((sign * fit.predict(list(x)), x) for x in product(*axes)),
                    reverse=True)
    best = None
    for _, start in scored[:5]:
        x = list(start)
        value = sign * fit.predict(x)
        for _ in range(200):
            improved = False
            for i in range(k):
                lo, hi = fit.bounds[i]

                def along(t, i=i):
                    z = list(x)
                    z[i] = t
                    return sign * fit.predict(z)
                cand = max(_cubic_on_interval(along, lo, hi), key=along)
                v = along(cand)
                if v > value + 1e-12 * max(1.0, abs(value)):
                    x[i], value, improved = cand, v, True
            if not improved:
                break
        if best is None or value > best[1]:
            best = (x, value)
    return best[0], sign * best[1]


def actual_equation(fit, centres, halves):
    """Coefficients in actual units, x = (X - c) / h: [(name, coefficient)] with names
    "Constant", "A", "A^2", "AB", "A^2B", ... (for an average block): every term is
    expanded as a product of (X_i - c_i) / h_i."""
    poly = {}
    for t in fit.terms:
        if t.kind == "block":
            continue
        expanded = {(): t.coef}
        for i in monomial(t):
            nxt = {}
            for m, v in expanded.items():
                for mm, vv in ((tuple(sorted(m + (i,))), v / halves[i]),
                               (m, -v * centres[i] / halves[i])):
                    nxt[mm] = nxt.get(mm, 0.0) + vv
            expanded = nxt
        for m, v in expanded.items():
            poly[m] = poly.get(m, 0.0) + v
    order = sorted(poly, key=lambda m: (len(m), len(set(m)) > 1, m))
    return [(monomial_name(m), poly[m]) for m in order]


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
