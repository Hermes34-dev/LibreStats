"""Per-observation regression diagnostics and Cook's distance (Cook, 1977).

For a full-rank least-squares fit y = X b + e with p parameters (intercept,
block and any other model columns included), n observations and
s^2 = SSE / (n - p):

- leverage           h_i = x_i' (X'X)^-1 x_i  (Cook's v_i)
- studentized        t_i = e_i / (s sqrt(1 - h_i))   (internally studentized)
- variance ratio     V(yhat_i) / V(e_i) = h_i / (1 - h_i)
- Cook's distance    D_i = (b_(i) - b)' X'X (b_(i) - b) / (p s^2)
                         = t_i^2 h_i / (p (1 - h_i))           (Cook's eq. 7)

D_i is read as a level of the F(p, n - p) distribution: removing observation i
moves the estimate to the edge of the 100 F_cdf(D_i)% confidence region for b
(Cook, 1977). An observation is flagged as influential when that level is 50%
or more, i.e. D_i >= the median of F(p, n - p) (Bollen & Jackman, 1990).

D_i is undefined when h_i = 1 (the fit without observation i is not full rank)
and t_i when s = 0 (a perfect fit); both are then None.

Outlier test (Cook & Weisberg, 1982; Weisberg, 2014): the externally studentized
residual (the residual scaled by s from the fit without observation i)
    t*_i = t_i sqrt((n - p - 1) / (n - p - t_i^2))
has a t distribution with n - p - 1 df if observation i is not an outlier;
p_i = 2 P(T > |t*_i|), and with a Bonferroni correction for testing all n
observations p_bonf = min(1, n p_i). Not available when n - p - 1 <= 0, or when
t_i^2 >= n - p (the other observations then fit exactly, t*_i is infinite).

Flags: "outlier" if p_bonf < 0.05; "check" if p_i < 0.05 (but not an outlier);
"high leverage" if h_i > 2p/n (Hoaglin & Welsch, 1978; not used when p/n > 0.5,
where every h_i could exceed it); "influential" as above.
"""

import math
from dataclasses import dataclass

from .distributions import f_cdf, t_sf

FLAG_LEVEL = 0.5          # flag D_i at or above the median of F(p, n - p)
OUTLIER_ALPHA = 0.05      # Bonferroni p below this: outlier; unadjusted p: check
LEVERAGE_ONE = 1e-10      # 1 - h_i below this: h_i = 1


@dataclass
class Diagnostic:
    observed: float
    fitted: float
    residual: float
    leverage: float
    studentized: float = None
    ratio: float = None        # V(yhat_i) / V(e_i) = h_i / (1 - h_i)
    cooks: float = None
    level: float = None        # F(p, n - p) cdf at D_i: the confidence-region level
    external: float = None     # externally studentized residual t*_i
    p_outlier: float = None    # two-sided p of t*_i, t(n - p - 1)
    p_bonferroni: float = None
    high_leverage: bool = False

    @property
    def influential(self):
        return self.level is not None and self.level >= FLAG_LEVEL

    @property
    def outlier(self):
        return self.p_bonferroni is not None and self.p_bonferroni < OUTLIER_ALPHA

    @property
    def check(self):
        return not self.outlier and self.p_outlier is not None and \
            self.p_outlier < OUTLIER_ALPHA

    def flags(self):
        return [name for name, on in (("outlier", self.outlier), ("check", self.check),
                                      ("high leverage", self.high_leverage),
                                      ("influential", self.influential)) if on]


def leverage_cutoff(n, p):
    """2p/n (Hoaglin & Welsch, 1978), or None when p/n > 0.5 (the rule cannot work)."""
    return None if p / n > 0.5 else 2 * p / n


def diagnostics(rows, y, beta, inv, sse):
    """Diagnostics of each row of the model matrix (rows), given the fit: beta, the
    inverse of X'X and the residual sum of squares. Needs n > p."""
    n, p = len(rows), len(rows[0])
    df = n - p
    if df <= 0:
        raise ValueError("Diagnostics need residual degrees of freedom (n > p).")
    s2 = sse / df
    cutoff = leverage_cutoff(n, p)
    out = []
    for r, v in zip(rows, y):
        fit = sum(a * b for a, b in zip(r, beta))
        e = v - fit
        h = sum(r[a] * sum(inv[a][c] * r[c] for c in range(p)) for a in range(p))
        h = min(max(h, 0.0), 1.0)
        d = Diagnostic(v, fit, e, h)
        if 1 - h > LEVERAGE_ONE:
            d.ratio = h / (1 - h)
            if s2 > 0:
                d.studentized = e / (s2 * (1 - h)) ** 0.5
                d.cooks = d.studentized ** 2 * d.ratio / p
                d.level = f_cdf(d.cooks, p, df)
                r2 = d.studentized ** 2
                if df - 1 > 0 and df - r2 > 1e-12 * df:
                    d.external = d.studentized * math.sqrt((df - 1) / (df - r2))
                    d.p_outlier = min(1.0, 2 * t_sf(abs(d.external), df - 1))
                    d.p_bonferroni = min(1.0, n * d.p_outlier)
        d.high_leverage = cutoff is not None and h > cutoff
        out.append(d)
    return out


def _g(x, digits=2):
    return f"{x:.{digits}g}" if abs(x) < 10 ** digits else f"{x:.0f}"


def verdict(label, d, n, p):
    """One line for a flagged observation: what was found and what to do; None if no flag."""
    if not d.flags():
        return None
    found, todo = [], []
    if d.outlier:
        found.append(f"outlier (t* = {d.external:.2f}, Bonferroni p = {d.p_bonferroni:.2g})")
        todo.append("check the record and how the run was done; if it was an error, correct "
                    "or exclude it and refit, otherwise keep it and report it")
    elif d.check:
        found.append(f"possible outlier (t* = {d.external:.2f}, p = {d.p_outlier:.2g}, but "
                     f"Bonferroni p = {d.p_bonferroni:.2g})")
        todo.append(f"check the record; among {n} runs about {0.05 * n:.2g} such p values "
                    "below 0.05 are expected by chance")
    if d.high_leverage:
        found.append(f"high leverage (h = {d.leverage:.3g} > 2p/n = {2 * p / n:.3g})")
        todo.append("the fit relies on this run's settings; consider replicating it")
    if d.influential:
        found.append(f"influential (D = {_g(d.cooks)}, F level {100 * d.level:.0f}%)")
        todo.append("compare the fit with and without this run before drawing conclusions")
    elif d.cooks is not None:
        found.append(f"not influential (D = {_g(d.cooks)}, F level {100 * d.level:.1f}%)")
    return f"Run {label}: " + "; ".join(found) + ". Action: " + "; ".join(todo) + "."
