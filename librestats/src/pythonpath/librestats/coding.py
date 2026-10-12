"""Coding and recoding: dummy variables, presence/absence, bins, recoding rules.

Categorical values may be numbers (float) or text; None or "" means missing.

References
----------
Wang, H., & Song, M. (2011). Ckmeans.1d.dp: optimal k-means clustering in one
    dimension by dynamic programming. The R Journal, 3(2), 29-33.
"""

import bisect
import math
from dataclasses import dataclass, field

from .outliers import quantile


class CodingError(ValueError):
    pass


def is_missing(v):
    return v is None or (isinstance(v, str) and not v.strip())


def _key(v):
    """Sort/compare key: numbers before text, numbers by value, text alphabetically."""
    return (0, v, "") if isinstance(v, float) else (1, 0.0, str(v).strip())


def label(v):
    """Display text for a category: 2.0 -> '2', 2.5 -> '2.5', text unchanged."""
    if isinstance(v, float):
        return str(int(v)) if v.is_integer() else repr(v)
    return str(v).strip()


def categories(values):
    """Distinct non-missing values in sorted order."""
    seen = {}
    for v in values:
        if not is_missing(v):
            seen.setdefault(_key(v), v)
    return [seen[k] for k in sorted(seen)]


@dataclass
class Coded:
    names: list                  # output column suffixes / names
    columns: list                # one list of output values per new column
    params: list = field(default_factory=list)
    table: list = field(default_factory=list)   # extra report rows (e.g. bins)
    notes: list = field(default_factory=list)


# --- Dummy variables and presence/absence ------------------------------------------

def dummies(values, drop_first=False):
    """One 0/1 column per category (one-hot); missing stays missing.

    With drop_first, the first category is the reference and gets no column
    (k - 1 dummies, as needed in regression with an intercept)."""
    cats = categories(values)
    if len(cats) < 2:
        raise CodingError("Dummy coding needs at least 2 categories.")
    used = cats[1:] if drop_first else cats
    columns = []
    for cat in used:
        k = _key(cat)
        columns.append([None if is_missing(v) else (1.0 if _key(v) == k else 0.0)
                        for v in values])
    params = [("Categories", len(cats))]
    if drop_first:
        params.append(("Reference (no column)", label(cats[0])))
    return Coded([label(c) for c in used], columns, params)


def presence_absence(values, threshold=0.0, empty_is_absent=False):
    """1 if a number is above the threshold, else 0. Text is treated as missing."""
    out, text = [], 0
    for v in values:
        if isinstance(v, float):
            out.append(1.0 if v > threshold else 0.0)
        elif is_missing(v):
            out.append(0.0 if empty_is_absent else None)
        else:
            out.append(None)
            text += 1
    params = [("Present if value >", threshold),
              ("Empty cells", "absent (0)" if empty_is_absent else "missing")]
    notes = [f"{text} text cell(s) left empty."] if text else []
    return Coded(["present"], [out], params, notes=notes)


# --- Discretisation -----------------------------------------------------------------

def _numbers(values):
    xs = [v for v in values if isinstance(v, float)]
    if len(xs) < 2:
        raise CodingError("Binning needs at least 2 numeric values.")
    return xs


def _fmt(x):
    return f"{x:.6g}"


def _codes_from_edges(values, edges):
    """Bin numbers 1..k for edges e0 < e1 < ... < ek; bins are (e(i-1), e(i)],
    the first bin also includes e0 (as pandas.cut / qcut with include_lowest)."""
    inner = edges[1:-1]
    return [float(bisect.bisect_left(inner, v) + 1) if isinstance(v, float) else None
            for v in values]


def _bin_table(values, codes, edges):
    rows = [["Bin", "From (>)", "To (<=)", "Count", "Label"]]
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        lab = f"{'[' if i == 0 else '('}{_fmt(lo)}, {_fmt(hi)}]"
        rows.append([i + 1, lo, hi, sum(1 for c in codes if c == i + 1), lab])
    return rows


def bins_equal_width(values, k):
    xs = _numbers(values)
    lo, hi = min(xs), max(xs)
    if lo == hi:
        raise CodingError("All values are equal.")
    step = (hi - lo) / k
    edges = [lo + i * step for i in range(k)] + [hi]
    codes = _codes_from_edges(values, edges)
    return Coded(["bin"], [codes], [("Method", "equal width"), ("Bins", k), ("Width", step)],
                 _bin_table(values, codes, edges))


def bins_quantile(values, k):
    xs = sorted(_numbers(values))
    edges = [quantile(xs, i / k) for i in range(k + 1)]
    unique = sorted(set(edges))
    notes = []
    if len(unique) < len(edges):
        notes.append(f"Tied values made some quantile boundaries equal: {len(unique) - 1} "
                     f"bins instead of {k}.")
    if len(unique) < 2:
        raise CodingError("All values are equal.")
    codes = _codes_from_edges(values, unique)
    return Coded(["bin"], [codes], [("Method", "quantiles (equal counts)"),
                                    ("Bins", len(unique) - 1)],
                 _bin_table(values, codes, unique), notes)


def kmeans_1d(xs, k):
    """Optimal 1-D k-means (minimum within-cluster sum of squares) by dynamic
    programming with divide-and-conquer optimisation (Wang & Song, 2011).
    Returns (sorted values, cluster index for each sorted value, total SS)."""
    x = sorted(xs)
    n = len(x)
    mean = sum(x) / n
    s1, s2 = [0.0], [0.0]
    for v in x:  # prefix sums of centred values for numerical stability
        s1.append(s1[-1] + (v - mean))
        s2.append(s2[-1] + (v - mean) ** 2)

    def cost(j, i):  # SS of x[j..i] inclusive
        m = i - j + 1
        a = s1[i + 1] - s1[j]
        return max(0.0, s2[i + 1] - s2[j] - a * a / m)

    prev = [cost(0, i) for i in range(n)]
    back = [[0] * n]
    for m in range(1, k):
        cur = [math.inf] * n
        arg = [0] * n

        def solve(lo, hi, jlo, jhi):
            if lo > hi:
                return
            mid = (lo + hi) // 2
            best, best_j = math.inf, jlo
            for j in range(max(jlo, m), min(jhi, mid) + 1):
                c = prev[j - 1] + cost(j, mid)
                if c < best:
                    best, best_j = c, j
            cur[mid], arg[mid] = best, best_j
            solve(lo, mid - 1, jlo, best_j)
            solve(mid + 1, hi, best_j, jhi)

        solve(m, n - 1, m, n - 1)
        prev = cur
        back.append(arg)
    labels = [0] * n
    i = n - 1
    for m in range(k - 1, -1, -1):
        j = back[m][i] if m > 0 else 0
        for t in range(j, i + 1):
            labels[t] = m
        i = j - 1
    return x, labels, prev[n - 1]


def bins_kmeans(values, k):
    xs = _numbers(values)
    distinct = len(set(xs))
    notes = []
    if distinct < k:
        notes.append(f"Only {distinct} distinct values: {distinct} bins instead of {k}.")
        k = distinct
    if k < 2:
        raise CodingError("All values are equal.")
    x, labels, ss = kmeans_1d(xs, k)
    upper = [max(v for v, lab in zip(x, labels) if lab == c) for c in range(k)]
    lower = [min(v for v, lab in zip(x, labels) if lab == c) for c in range(k)]
    codes = [None if not isinstance(v, float) else float(bisect.bisect_left(upper, v) + 1)
             for v in values]
    table = [["Bin", "Smallest", "Largest", "Count", "Centre"]]
    for c in range(k):
        members = [v for v, lab in zip(x, labels) if lab == c]
        table.append([c + 1, lower[c], upper[c], len(members), sum(members) / len(members)])
    total = sum((v - sum(x) / len(x)) ** 2 for v in x)
    return Coded(["bin"], [codes], [("Method", "k-means (optimal, 1-D)"), ("Bins", k),
                                    ("Between-bin SS / total SS",
                                     1 - ss / total if total else 1.0)],
                 table, notes)


# --- Recoding -----------------------------------------------------------------------

@dataclass
class Rules:
    mapping: list            # [(set of keys, new value)]
    default: object = None   # value for unmatched categories, or None to keep them


def _parse_value(text):
    text = text.strip()
    try:
        return float(text.replace(",", ".")) if text else ""
    except ValueError:
        return text


def parse_rules(text):
    """Parse lines like  '1; 2 -> low'  or  'yes; y -> 1'  or  '* -> other'.

    Several old values are separated by ';'. '*' matches every value not
    listed elsewhere. Blank lines and lines starting with '#' are ignored."""
    mapping, default = [], None
    for number, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "->" not in line:
            raise CodingError(f"Line {number}: expected 'old -> new' (got {line!r}).")
        left, right = line.rsplit("->", 1)
        new = _parse_value(right)
        olds = [o.strip() for o in left.split(";") if o.strip()]
        if not olds:
            raise CodingError(f"Line {number}: no old value before '->'.")
        if olds == ["*"]:
            default = new
            continue
        mapping.append(({_key(_parse_value(o)) for o in olds}, new))
    if not mapping and default is None:
        raise CodingError("No recoding rules given.")
    return Rules(mapping, default)


def recode(values, rules):
    out, unmatched = [], set()
    for v in values:
        if is_missing(v):
            out.append(None)
            continue
        k = _key(v)
        for keys, new in rules.mapping:
            if k in keys:
                out.append(new)
                break
        else:
            if rules.default is not None:
                out.append(rules.default)
            else:
                out.append(v)
                unmatched.add(label(v))
    notes = []
    if unmatched:
        notes.append("Kept unchanged (no rule): " + ", ".join(sorted(unmatched)))
    return Coded(["recoded"], [out], [("Rules", len(rules.mapping))], notes=notes)


def rules_template(values):
    """'value -> value' lines for every category, as a starting point for editing."""
    return "\n".join(f"{label(c)} -> {label(c)}" for c in categories(values))
