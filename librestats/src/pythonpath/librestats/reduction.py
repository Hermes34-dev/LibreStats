"""Model reduction for response surface models: manual, backward elimination by
p-value, and all subsets ranked by AICc.

Terms are monomials of the coded factors, e.g. A^2B = A.A.B. The parents of a
term are the lower-order terms that divide it (A^2B -> A, B, A^2, AB). Only
parents that are in the starting model count: in a reduced cubic some parents
may have been removed as aliases. The constant and block terms are never removed.

- Hierarchical model: every term's parents (in the starting model) are in it.
- Removable term (hierarchy kept): no other term in the model has it as a parent.
- Locked terms are never removed; with hierarchy kept, neither are their parents.

Backward elimination: among the removable, unlocked terms take the one with the
highest p; if p > alpha-out remove it, refit and repeat (one term per step, since
the p-values change after every refit).

All subsets: every subset of the candidate terms (hierarchical when hierarchy is
kept, always with the locked terms) is fitted and ranked by
    AICc = n ln(SSE/n) + 2p + 2p(p+1)/(n-p-1)      (p includes the constant)
lowest first; ties by fewer terms, then by higher predicted R^2. AICc only needs
SSE, which comes from the cross-product matrix of the starting model; R^2,
predicted R^2 and lack of fit are computed for the models shown.
"""

import math
from dataclasses import dataclass, field
from itertools import combinations

from . import rsm
from .linalg import SingularMatrixError, cholesky, cho_solve

METHODS = ("none", "manual", "backward", "subsets")
METHOD_NAMES = {"none": "none", "manual": "manual", "backward": "backward elimination",
                "subsets": "all subsets (AICc)"}
MAX_CANDIDATES = 20       # all subsets: more candidate terms -> use backward elimination
MAX_SUBSETS = 200000      # all subsets: more models to fit -> refused (about 10 s)
PRED_LEVERAGE = 0.9999    # all-subsets table: predicted R^2 n/a at or above this leverage


class ReductionError(ValueError):
    pass


@dataclass
class Candidate:
    """One model of the all-subsets search."""
    terms: list               # term names (without the constant and blocks)
    p: int                    # parameters including the constant and blocks
    aicc: float
    r2: float = None
    r2_adj: float = None
    r2_pred: float = None
    lof_p: float = None


@dataclass
class Reduction:
    method: str
    alpha: float = None
    hierarchy: bool = True
    locked: list = field(default_factory=list)
    start: object = None              # rsm.Fit of the starting model
    final: object = None              # rsm.Fit of the final model
    removed: list = field(default_factory=list)       # term names, in order
    steps: list = field(default_factory=list)         # backward: (term, p)
    candidates: list = field(default_factory=list)    # all subsets: top N (Candidate)
    chosen: object = None                             # all subsets: the Candidate used
    ranked: list = field(default_factory=list)        # all subsets: ranked so far
    evaluated: int = 0                                # all subsets: models fitted
    rank: int = 1                                     # all subsets: the model used
    warnings: list = field(default_factory=list)
    excluded: list = field(default_factory=list)      # aliased terms left out of the start

    def summary(self, model_name):
        how = {"manual": "manually", "backward": "by backward elimination",
               "subsets": "by all subsets (AICc)"}[self.method]
        opts = []
        if self.method == "backward":
            opts.append(f"alpha = {self.alpha:g}")
        if self.method == "subsets":
            opts.append(f"rank {self.rank} of {self.evaluated} models")
        opts.append("hierarchy kept" if self.hierarchy else "hierarchy off")
        if self.locked:
            opts.append("locked: " + ", ".join(self.locked))
        removed = ", ".join(self.removed) if self.removed else "none"
        return (f"Reduced from {model_name} {how} ({'; '.join(opts)}). "
                f"Removed: {removed}.")


def parents(term_name):
    """Names of the lower-order terms that divide a term: 'A^2B' -> A, B, A^2, AB."""
    m = rsm.name_monomial(term_name)
    out = set()
    for size in range(1, len(m)):
        for sub in combinations(m, size):
            out.add(rsm.monomial_name(tuple(sorted(sub))))
    return out


def missing_parents(names, start):
    """{term: parents in the starting model but not in names}, for terms in names."""
    names, start = set(names), set(start)
    out = {}
    for t in names:
        gone = (parents(t) & start) - names
        if gone:
            out[t] = sorted(gone, key=_order)
    return out


def removable(names, start, locked, hierarchy=True):
    """Terms of the model (names) that may be removed."""
    names = set(names)
    out = []
    for t in names:
        if t in locked:
            continue
        if hierarchy and any(t in parents(u) for u in names if u != t):
            continue
        out.append(t)
    return sorted(out, key=_order)


def _order(name):
    m = rsm.name_monomial(name)
    return (len(m), len(set(m)) > 1, m)


def _model_terms(fit):
    return [t.name for t in fit.terms if t.kind not in ("constant", "block")]


def reduce(coded, blocks, y, model, method="none", alpha=0.05, hierarchy=True, locked=(),
           remove=(), top=5, rank=1, exclude_aliased=False):
    """Fit the starting model and reduce it. Returns a Reduction whose .final is the
    rsm.Fit of the final model. exclude_aliased: leave the terms that stand in for an
    aliased term (rsm.alias_marks, e.g. AB^2 aliased with AC^2) out of the starting
    model, so the selection cannot pick a term whose meaning is ambiguous."""
    if method not in METHODS:
        raise ReductionError(f"Unknown model reduction {method!r}.")
    start = rsm.fit_model(coded, blocks, y, model)
    excluded = []
    if exclude_aliased and start.aliases:
        marks = rsm.alias_marks(start.aliases)
        excluded = sorted((t for t in _model_terms(start) if t in marks), key=_order)
        locked_marked = [t for t in locked if t.strip() in marks]
        if locked_marked:
            raise ReductionError(f"{locked_marked[0].strip()} is both locked and excluded as "
                                 "an aliased term.")
        if excluded:
            full_aliases = start.aliases
            start = rsm.fit_model(coded, blocks, y, model,
                                  keep=[t for t in _model_terms(start) if t not in marks])
            start.aliases = full_aliases
    names = _model_terms(start)
    start_names = list(names)
    aliased = {a.term for a in start.aliases}
    locked = [t.strip() for t in locked if t.strip()]
    for t in locked:
        if t in aliased:
            raise ReductionError(f"{t} cannot be locked: it is aliased with other terms in "
                                 "this design and was removed from the reduced cubic.")
        if t not in names:
            raise ReductionError(f"{t} cannot be locked: it is not a term of the model.")
    red = Reduction(method, alpha if method == "backward" else None, hierarchy, locked,
                    start)
    red.excluded = excluded
    if hierarchy:
        # parents of locked terms cannot be removed either
        locked_all = set(locked)
        for t in locked:
            locked_all |= parents(t) & set(names)
    else:
        locked_all = set(locked)

    if method == "none":
        red.final = start
    elif method == "manual":
        remove = [t.strip() for t in remove if t.strip()]
        for t in remove:
            if t in locked:
                raise ReductionError(f"{t} is locked and cannot be removed.")
            if t not in names and t not in aliased:
                raise ReductionError(f"{t} is not a term of the model.")
        gone = [t for t in remove if t in names]
        keep = [t for t in names if t not in gone]
        if hierarchy:
            for t in sorted(gone, key=_order):
                children = sorted((u for u in keep if t in parents(u)), key=_order)
                if children:
                    listed = children[0] if len(children) == 1 else \
                        ", ".join(children[:-1]) + " and " + children[-1]
                    verb = "is" if len(children) == 1 else "are"
                    raise ReductionError(
                        f"Keeping the model hierarchical: {t} cannot be removed while "
                        f"{listed} {verb} in the model. Remove "
                        f"{'it' if len(children) == 1 else 'them'} too, or untick 'Keep "
                        "model hierarchy'.")
        red.removed = sorted(gone, key=_order)
        red.final = rsm.fit_model(coded, blocks, y, model, keep=keep)
        if [t for t in remove if t in aliased]:
            red.warnings.append("Aliased terms are removed from the reduced cubic anyway: "
                                + ", ".join(t for t in remove if t in aliased) + ".")
    elif method == "backward":
        current, fit = list(names), start
        while True:
            options = removable(current, start_names, locked_all, hierarchy)
            pvals = {t.name: t.p for t in fit.terms if t.name in options}
            if not pvals:
                break
            worst = max(pvals, key=lambda t: (pvals[t], _order(t)))
            if pvals[worst] <= alpha:
                break
            red.steps.append((worst, pvals[worst]))
            red.removed.append(worst)
            current.remove(worst)
            fit = rsm.fit_model(coded, blocks, y, model, keep=current)
        red.final = fit
    else:
        _subsets(red, coded, blocks, y, model, names, start_names, locked_all, hierarchy,
                 top, rank)
    if not hierarchy:
        broken = missing_parents(_model_terms(red.final), start_names)
        text = ("Model hierarchy is off: the final model may contain a term without its "
                "lower-order terms, which makes the model depend on how the factors are "
                "coded.")
        if broken:
            text += " Missing: " + "; ".join(f"{t} without {', '.join(ps)}"
                                             for t, ps in sorted(broken.items())) + "."
        red.warnings.append(text)
    return red


def _subsets(red, coded, blocks, y, model, names, start_names, locked, hierarchy, top,
             rank):
    candidates = sorted((t for t in names if t not in locked), key=_order)
    if len(candidates) > MAX_CANDIDATES:
        raise ReductionError(f"All subsets: {len(candidates)} candidate terms is too many "
                             f"(at most {MAX_CANDIDATES}); use backward elimination, or "
                             "lock terms.")
    fixed = [t for t in names if t in locked]
    start = red.start
    # cross products of the starting model (constant, blocks and all terms)
    obs = start.observations
    mean = sum(y[i] for i in obs) / len(obs)
    yy = [y[i] - mean for i in obs]     # exact: the constant is always in the model
    cols = [t for t in start.terms]
    rows = [[rsm._value(t, coded[i]) if t.kind != "block" else _block_value(t, blocks, i,
                                                                             obs)
             for t in cols] for i in obs]
    p_all = len(cols)
    xtx = [[sum(r[a] * r[b] for r in rows) for b in range(p_all)] for a in range(p_all)]
    xty = [sum(r[a] * v for r, v in zip(rows, yy)) for a in range(p_all)]
    yty = sum(v * v for v in yy)
    n = len(yy)
    index = {t.name: j for j, t in enumerate(cols)}
    base = [j for j, t in enumerate(cols) if t.kind in ("constant", "block")] + \
        [index[t] for t in fixed]

    found = []

    def evaluate(chosen):
        idx = base + [index[t] for t in chosen]
        p = len(idx)
        if n - p - 1 <= 0:
            return
        try:
            L = cholesky([[xtx[a][b] for b in idx] for a in idx])
        except SingularMatrixError:
            return
        b = cho_solve(L, [xty[a] for a in idx])
        sse = max(yty - sum(bi * xty[a] for bi, a in zip(b, idx)), 0.0)
        value = rsm.aicc(sse, n, p)
        if value is not None:
            found.append(Candidate(list(chosen), p, value))

    count = [0]

    def walk(i, chosen):
        if i == len(candidates):
            count[0] += 1
            if count[0] > MAX_SUBSETS:
                raise ReductionError(
                    f"All subsets: more than {MAX_SUBSETS} models to fit. Keep the model "
                    "hierarchy, lock terms, or use backward elimination.")
            evaluate(chosen)
            return
        t = candidates[i]
        walk(i + 1, chosen)                                     # without t
        have = set(chosen) | set(fixed)
        if not hierarchy or (parents(t) & set(start_names)) <= have:
            walk(i + 1, chosen + [t])                           # with t

    walk(0, [])
    red.evaluated = len(found)
    if not found:
        raise ReductionError("All subsets: no model leaves residual degrees of freedom "
                             "(n - p - 1 > 0).")
    found.sort(key=lambda c: (round(c.aicc, 9), c.p))
    # predicted R^2 breaks remaining ties; computed for the models near the top
    shown = []
    i = 0
    while i < len(found) and len(shown) < max(top, rank):
        j = i
        while j < len(found) and (round(found[j].aicc, 9), found[j].p) == \
                (round(found[i].aicc, 9), found[i].p):
            j += 1
        group = found[i:j]
        for c in group:
            _complete(c, coded, blocks, y, model, fixed)
        group.sort(key=lambda c: -(c.r2_pred if c.r2_pred is not None else -math.inf))
        shown += group
        i = j
    red.candidates = shown[:top]
    red.ranked = shown
    if not 1 <= rank <= len(shown):
        raise ReductionError(f"All subsets: rank {rank} is not available "
                             f"({len(shown)} models).")
    red.rank = rank
    chosen = red.chosen = shown[rank - 1]
    keep = fixed + chosen.terms
    red.removed = [t for t in sorted(names, key=_order) if t not in keep]
    red.final = rsm.fit_model(coded, blocks, y, model, keep=keep)


def _block_value(term, blocks, i, obs):
    levels = sorted({blocks[j] for j in obs})
    b = blocks[i]
    return 1.0 if b == levels[term.factors[0]] else -1.0 if b == levels[-1] else 0.0


def _complete(c, coded, blocks, y, model, fixed):
    fit = rsm.fit_model(coded, blocks, y, model, keep=fixed + c.terms)
    c.r2, c.r2_adj, c.lof_p = fit.r2, fit.r2_adj, fit.lof_p
    c.r2_pred = fit.r2_pred
    if any(d.leverage >= PRED_LEVERAGE for d in fit.diagnostics):
        c.r2_pred = None
