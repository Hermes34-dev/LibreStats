"""Analysis of two-level screening designs.

The factorial runs are coded -1/+1, centre points 0. Effects are estimated by
least squares, so a missing response does not stop the analysis (it does make
the effects slightly correlated). The model holds an intercept, block terms, a
centre-point term (when there are centre points, so that curvature does not bias
the intercept or the error) and the main effects plus, by default, the
two-factor interactions. An effect is twice its coefficient: the change in the
mean response from the low to the high level.

Aliasing is read from the data: two effect columns that are equal (or opposite)
over the factorial runs are aliased and estimated once; columns that are
correlated but not equal are partially aliased (non-regular designs such as
Plackett-Burman) and are left out of the model.

- Error: the residual of the model; if the model is saturated (no residual
  degrees of freedom), Lenth's pseudo standard error.
- Lenth (1989): PSE = 1.5 median{|c| : |c| < 2.5 s0}, s0 = 1.5 median |c|,
  over all m estimable contrasts; margin of error ME = t(1 - alpha/2; m/3) PSE,
  simultaneous margin SME = t(gamma; m/3) PSE, gamma = (1 + (1 - alpha)^(1/m))/2.
- Half-normal plot (Daniel, 1959): sorted |effects| against
  Phi^-1(0.5 + 0.5 (i - 0.5) / m).
- Curvature (Montgomery, 2017, section 6.8): the centre-point term estimates
  mean(centre) - mean(factorial); SS = estimate^2 / its variance factor (equal to
  nF nC (yF - yC)^2 / (nF + nC) without blocks), F test against pure error.
"""

import math
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations

from .distributions import f_sf, norm_ppf, t_isf, t_sf
from .doe import LETTERS
from .linalg import SingularMatrixError, inverse_spd, matvec

MODELS = ("2fi", "main")
RSM_POINT_TYPES = (-1, 2)  # axial points (central composite), edge midpoints (Box-Behnken)
RSM_DESIGN = ("This is a response surface design (central composite or Box-Behnken: it has "
              "axial or edge points). Analyse Screening Design works only on two-level "
              "screening designs (full or fractional factorial, Plackett-Burman). Use DOE / "
              "RSM > Analyse Response Surface Design instead.")
SHOW_ORDER = 3           # alias members up to this order are named


class ScreeningError(ValueError):
    pass


@dataclass
class Term:
    mask: int                # first (lowest-order) member of the alias chain
    aliases: list            # (sign, mask) of the other members, in order
    effect: float = None
    coef: float = None
    se: float = None         # standard error of the effect
    t: float = None
    p: float = None

    @property
    def name(self):
        return effect_name(self.mask)

    def alias_text(self, max_order=SHOW_ORDER):
        shown = [(s, m) for s, m in self.aliases if _order(m) <= max_order]
        more = len(self.aliases) - len(shown)
        text = " ".join(("+ " if s > 0 else "- ") + effect_name(m) for s, m in shown)
        if more:
            text += f" (+ {more} higher-order)"
        return text.lstrip("+ ") if text.startswith("+") else text


@dataclass
class Lenth:
    m: int
    pse: float
    me: float
    sme: float
    df: float


@dataclass
class Curvature:
    mean_factorial: float
    mean_centre: float
    n_factorial: int
    n_centre: int
    estimate: float          # centre minus factorial (block-adjusted)
    ss: float
    f: float = None
    p: float = None


@dataclass
class Analysis:
    n: int                               # runs with a response
    missing: int
    terms: list                          # Term, in model order
    intercept: float
    df: int                              # residual degrees of freedom
    mse: float
    r2: float
    r2_adj: float
    error: str                           # "residual" | "lenth" | "none"
    pure_df: int
    pure_ss: float
    contrasts: list                      # (Term, effect) of all estimable contrasts
    lenth: Lenth = None
    curvature: Curvature = None
    block_effects: list = field(default_factory=list)    # chains confounded with blocks
    notes: list = field(default_factory=list)

    def critical_t(self, alpha):
        if self.error == "residual":
            return t_isf(alpha / 2, self.df)
        if self.error == "lenth":
            return t_isf(alpha / 2, self.lenth.df)
        return None


def effect_name(mask):
    return "".join(LETTERS[i] for i in range(len(LETTERS)) if mask >> i & 1)


def _order(mask):
    return bin(mask).count("1")


def _median(values):
    v = sorted(values)
    n = len(v)
    return (v[n // 2] + v[(n - 1) // 2]) / 2


def lenth(effects, alpha=0.05):
    """Lenth's (1989) pseudo standard error and margins of error."""
    m = len(effects)
    a = [abs(e) for e in effects]
    s0 = 1.5 * _median(a)
    trimmed = [v for v in a if v < 2.5 * s0]
    pse = 1.5 * _median(trimmed) if trimmed else 0.0
    d = m / 3
    gamma = (1 + (1 - alpha) ** (1 / m)) / 2
    return Lenth(m, pse, t_isf(alpha / 2, d) * pse, t_isf(1 - gamma, d) * pse, d)


def half_normal_quantiles(m):
    return [norm_ppf(0.5 + 0.5 * (i - 0.5) / m) for i in range(1, m + 1)]


def _chains(coded, rows, blocks):
    """Alias chains of the effects over the factorial runs: (estimable chains,
    chains confounded with blocks, number of partially aliased two-factor interactions).
    All interactions are visited (up to order 2 for more than 12 factors), so every
    chain lists its aliases."""
    k = len(coded[0])
    block_rows = defaultdict(list)
    for pos, i in enumerate(rows):
        block_rows[blocks[i]].append(pos)
    blocked = len(block_rows) > 1
    max_order = k if (1 << k) <= 4096 else 2
    found = {}                       # column -> chain index
    chains, block_chain = [], []     # [(mask, aliases)], bool per chain
    partial = 0
    for order in range(1, max_order + 1):
        for members in combinations(range(k), order):
            mask = sum(1 << i for i in members)
            col = tuple(math.prod(coded[i][j] for j in members) for i in rows)
            neg = tuple(-v for v in col)
            if col in found or neg in found:
                idx = found[col] if col in found else found[neg]
                chains[idx][1].append((1 if col in found else -1, mask))
                continue
            if blocked and all(len({col[p] for p in pos}) == 1 for pos in block_rows.values()):
                found[col] = len(chains)
                chains.append((mask, []))
                block_chain.append(True)
                continue
            orthogonal = all(sum(col[p] for p in pos) == 0 for pos in block_rows.values()) \
                and all(sum(a * b for a, b in zip(col, rep)) == 0
                        for rep, idx in found.items() if not block_chain[idx])
            if orthogonal:
                found[col] = len(chains)
                chains.append((mask, []))
                block_chain.append(False)
            elif order == 2:
                partial += 1
    estimable = [Term(m, a) for (m, a), b in zip(chains, block_chain) if not b]
    confounded = [Term(m, a) for (m, a), b in zip(chains, block_chain) if b]
    return estimable, confounded, partial


def analyse(coded, ptype, blocks, y, model="2fi", alpha=0.05):
    """coded: rows of -1/0/+1; ptype: 1 factorial, 0 centre; y: responses (None = missing)."""
    if model not in MODELS:
        raise ScreeningError(f"Unknown model {model!r}.")
    if any(p in RSM_POINT_TYPES for p in ptype):
        raise ScreeningError(RSM_DESIGN)
    if any(p not in (0, 1) for p in ptype):
        raise ScreeningError("Point type must be 1 (factorial run) or 0 (centre point).")
    factorial = [i for i, p in enumerate(ptype) if p == 1]
    if any(v not in (-1, 1) for i in factorial for v in coded[i]):
        raise ScreeningError("The coded factor columns of the factorial runs must be -1 or +1.")
    if any(v != 0 for i, p in enumerate(ptype) if p == 0 for v in coded[i]):
        raise ScreeningError("Centre points must be coded 0.")
    obs = [i for i in range(len(y)) if y[i] is not None]
    missing = len(y) - len(obs)
    if len([i for i in obs if ptype[i] == 1]) < 3:
        raise ScreeningError("Enter the responses: at least 3 factorial runs need a value.")

    chains, confounded, partial = _chains(coded, factorial, blocks)
    notes = []
    if partial:
        notes.append(f"{partial} two-factor interaction(s) are partially aliased with other "
                     "effects (a non-regular design such as Plackett-Burman); they are left "
                     "out of the model.")
    if missing:
        notes.append(f"{missing} run(s) without a response are left out; the effect "
                     "estimates are then slightly correlated.")
    max_order = 2 if model == "2fi" else 1
    terms = [Term(t.mask, t.aliases) for t in chains if _order(t.mask) <= max_order]

    # least squares: intercept, blocks, centre-point term, effects
    obs_blocks = sorted({blocks[i] for i in obs})
    has_centre = any(ptype[i] == 0 for i in obs) and any(ptype[i] == 1 for i in obs)
    rows = []
    for i in obs:
        row = [1.0] + [1.0 if blocks[i] == b else 0.0 for b in obs_blocks[1:]]
        if has_centre:
            row.append(1.0 if ptype[i] == 0 else 0.0)
        row += [float(math.prod(coded[i][j] for j in range(len(coded[i])) if t.mask >> j & 1))
                for t in terms]
        rows.append(row)
    yy = [y[i] for i in obs]
    n, p = len(rows), len(rows[0])
    if p > n:
        raise ScreeningError(f"{n} runs with a response cannot estimate {p} model terms. "
                             "Choose main effects only.")
    xtx = [[sum(r[a] * r[b] for r in rows) for b in range(p)] for a in range(p)]
    xty = [sum(r[a] * v for r, v in zip(rows, yy)) for a in range(p)]
    try:
        inv = inverse_spd(xtx)
    except SingularMatrixError:
        raise ScreeningError("The runs with a response cannot separate all model terms "
                             "(too many missing responses). Choose main effects only.") from None
    beta = matvec(inv, xty)
    fitted = [sum(a * b for a, b in zip(r, beta)) for r in rows]
    sse = sum((v - f) ** 2 for v, f in zip(yy, fitted))
    df = n - p
    mean = sum(yy) / n
    sst = sum((v - mean) ** 2 for v in yy)
    mse = sse / df if df else None
    r2 = 1 - sse / sst if sst > 0 else None
    r2_adj = 1 - mse / (sst / (n - 1)) if sst > 0 and df else None

    # pure error: repeated points within a block
    groups = defaultdict(list)
    for i in obs:
        groups[(blocks[i], tuple(coded[i]))].append(y[i])
    pure_ss = sum(sum((v - sum(g) / len(g)) ** 2 for v in g) for g in groups.values())
    pure_df = sum(len(g) - 1 for g in groups.values())

    # all estimable contrasts (complete factorial data): effects for Lenth / half-normal
    first = 1 + len(obs_blocks) - 1 + (1 if has_centre else 0)
    contrasts, lenth_result = [], None
    complete = all(y[i] is not None for i in factorial)
    if complete:
        nf = len(factorial)
        for t in chains:
            col = [math.prod(coded[i][j] for j in range(len(coded[i])) if t.mask >> j & 1)
                   for i in factorial]
            contrasts.append((t, 2 * sum(c * y[i] for c, i in zip(col, factorial)) / nf))
        if len(contrasts) >= 3:
            lenth_result = lenth([e for _, e in contrasts], alpha)
    elif df == 0:
        notes.append("Lenth's method needs a response for every factorial run.")

    error = "residual" if df > 0 else "lenth" if lenth_result and lenth_result.pse > 0 \
        else "none"
    for j, t in enumerate(terms):
        b = beta[first + j]
        t.coef, t.effect = b, 2 * b
        if error == "residual":
            t.se = 2 * math.sqrt(mse * inv[first + j][first + j])
            t.t = t.effect / t.se if t.se > 0 else math.inf
            t.p = 2 * t_sf(abs(t.t), df)
        elif error == "lenth":
            t.se = lenth_result.pse
            t.t = t.effect / t.se
            t.p = 2 * t_sf(abs(t.t), lenth_result.df)
    if error == "none":
        notes.append("No error estimate: the model is saturated and Lenth's method is not "
                     "available; no t or p values.")

    curvature = None
    if has_centre:
        c = first - 1
        f_rows = [y[i] for i in obs if ptype[i] == 1]
        c_rows = [y[i] for i in obs if ptype[i] == 0]
        est = beta[c]
        ss = est ** 2 / inv[c][c]
        curvature = Curvature(sum(f_rows) / len(f_rows), sum(c_rows) / len(c_rows),
                              len(f_rows), len(c_rows), est, ss)
        if pure_df:
            curvature.f = ss / (pure_ss / pure_df)
            curvature.p = f_sf(curvature.f, 1, pure_df)
        else:
            notes.append("Curvature cannot be tested: it needs at least 2 centre points (or "
                         "replicated runs) for pure error.")

    return Analysis(n, missing, terms, beta[0], df, mse, r2, r2_adj, error, pure_df,
                    pure_ss, contrasts, lenth_result, curvature, confounded, notes)
