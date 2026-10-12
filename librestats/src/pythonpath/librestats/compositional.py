"""Log-ratio transforms for compositional data (pure Python).

Rows are compositions (e.g. fatty-acid percentages of one sample); columns are
parts. Log-ratios are scale invariant, so rows need not be closed to exactly
100 %, but every part must be positive.

References
----------
Aitchison, J. (1982). The statistical analysis of compositional data. JRSS B, 44(2), 139-177.
Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barcelo-Vidal, C. (2003).
    Isometric logratio transformations for compositional data analysis.
    Math. Geol., 35(3), 279-300.
Filzmoser, P., Hron, K., & Templ, M. (2018). Applied Compositional Data Analysis.
    Springer.  (pivot coordinates)
Martin-Fernandez, J. A., Barcelo-Vidal, C., & Pawlowsky-Glahn, V. (2003). Dealing with
    zeros and missing values in compositional data sets using nonparametric
    imputation. Math. Geol., 35(3), 253-278.  (multiplicative replacement)
"""

import math
from dataclasses import dataclass, field


class CompositionError(ValueError):
    pass


@dataclass
class LogRatioResult:
    names: list                 # output column names
    rows: list                  # one list per input row, or None for skipped rows
    params: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def _gmean_log(xs):
    return sum(math.log(x) for x in xs) / len(xs)


def clr(x):
    """Centred log-ratio: ln(x_i / g(x)), g = geometric mean of all parts."""
    g = _gmean_log(x)
    return [math.log(v) - g for v in x]


def ilr_pivot(x):
    """Isometric log-ratio, pivot coordinates (Filzmoser et al., 2018):
    z_i = sqrt((D-i)/(D-i+1)) * ln(x_i / g(x_{i+1}, ..., x_D)),  i = 1..D-1.
    z_1 holds all relative information about part 1 versus the other parts."""
    D = len(x)
    out = []
    for i in range(D - 1):
        rest = D - i - 1
        out.append(math.sqrt(rest / (rest + 1)) * (math.log(x[i]) - _gmean_log(x[i + 1:])))
    return out


def pivot_basis(D):
    """Orthonormal clr-space basis (D-1 rows) with ilr_pivot(x) = basis . clr(x)."""
    rows = []
    for i in range(D - 1):
        rest = D - i - 1
        f = math.sqrt(rest / (rest + 1))
        rows.append([0.0] * i + [f] + [-f / rest] * rest)
    return rows


def replace_zeros(rows, factor=0.65):
    """Multiplicative replacement of zeros (Martin-Fernandez et al., 2003).

    A zero in part j becomes delta_j = factor x (smallest positive value of part j),
    a proxy for the detection limit; the non-zero parts of that row are shrunk
    so the row total is unchanged and their ratios are preserved.
    Returns (new rows, deltas, number of zeros replaced).
    """
    D = len(rows[0])
    deltas = []
    for j in range(D):
        positive = [r[j] for r in rows if r is not None and r[j] > 0]
        if not positive:
            raise CompositionError(f"Part {j + 1} has no positive values.")
        deltas.append(factor * min(positive))
    out, replaced = [], 0
    for r in rows:
        if r is None:
            out.append(None)
            continue
        zeros = [j for j, v in enumerate(r) if v == 0]
        if not zeros:
            out.append(list(r))
            continue
        total = sum(r)
        added = sum(deltas[j] for j in zeros)
        if added >= total:
            raise CompositionError("A row has so many zeros that replacement is impossible.")
        shrink = 1 - added / total
        out.append([deltas[j] if v == 0 else v * shrink for j, v in enumerate(r)])
        replaced += len(zeros)
    return out, deltas, replaced


def log_ratios(rows, names, method="clr", zeros="error", factor=0.65):
    """CLR or ILR (pivot) coordinates for rows of parts (None = missing).

    Rows with a missing part are skipped. zeros: "error" or "replace".
    """
    D = len(names)
    if D < 2:
        raise CompositionError("A composition needs at least 2 parts.")
    clean = [None if any(v is None for v in r) else list(r) for r in rows]
    skipped = sum(1 for r in clean if r is None)
    if all(r is None for r in clean):
        raise CompositionError("No complete rows.")
    negative = sum(1 for r in clean if r is not None for v in r if v < 0)
    if negative:
        raise CompositionError(f"{negative} negative value(s): compositions must be non-negative.")
    params, notes = [("Parts (D)", D)], []
    n_zero = sum(1 for r in clean if r is not None for v in r if v == 0)
    if n_zero:
        if zeros != "replace":
            raise CompositionError(f"{n_zero} zero value(s): log-ratios need positive parts. "
                                   "Choose zero replacement, or remove the zeros.")
        clean, _, replaced = replace_zeros(clean, factor)
        params.append(("Zeros replaced", replaced))
        notes.append(f"Zeros replaced by {factor:g} x the smallest positive value of their part "
                     "(multiplicative replacement, Martin-Fernandez et al., 2003).")

    sums = [sum(r) for r in clean if r is not None]
    params += [("Smallest row total", min(sums)), ("Largest row total", max(sums))]
    if max(sums) - min(sums) > 1e-6 * max(sums):
        notes.append("Row totals differ. Log-ratios do not depend on the total, so the data "
                     "need not be closed; check that each row really is one composition.")
    if skipped:
        notes.append(f"{skipped} row(s) with a missing part were skipped.")

    if method == "clr":
        out_names = [f"clr({n})" for n in names]
        out = [None if r is None else clr(r) for r in clean]
        notes.append("CLR values sum to 0 in each row, so their covariance matrix is singular; "
                     "use ILR coordinates for methods that need a full-rank covariance.")
    elif method == "ilr":
        out_names = [f"ilr{i + 1}: {names[i]} vs ({', '.join(names[i + 1:])})"
                     for i in range(D - 1)]
        out = [None if r is None else ilr_pivot(r) for r in clean]
        notes.append("Pivot coordinates: ilr k compares part k with the geometric mean of the "
                     "parts after it. Reorder the columns to put the part of interest first.")
    else:
        raise ValueError(f"Unknown method {method!r}")
    params.insert(0, ("Method", "centred log-ratio (CLR)" if method == "clr"
                      else "isometric log-ratio (ILR, pivot coordinates)"))
    return LogRatioResult(out_names, out, params, notes)
