"""Histogram binning: frequency tables for numeric data.

Bins are half-open intervals [a, b); the last bin also includes its upper
limit, so the maximum is counted (as numpy.histogram). The number of bins is
chosen automatically with Sturges' rule (Sturges, 1926), k = ceil(log2 n) + 1,
or set by the user as a number of bins or a bin width.
"""

import bisect
import math
from dataclasses import dataclass, field

MAX_BINS = 1000


class HistogramError(ValueError):
    pass


def numeric_values(data):
    """Flatten a 2-D data array and keep the numbers; also count non-empty text cells."""
    values, text = [], 0
    for row in data:
        for v in row:
            if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v):
                values.append(float(v))
            elif isinstance(v, str) and v.strip():
                text += 1
    return values, text


def fmt(x):
    """Edge value for a bin label: up to 10 significant digits, no '-0'."""
    return f"{x:.10g}" if x != 0 else "0"


def _clean(x):
    """Remove floating-point noise such as 0.30000000000000004 -> 0.3."""
    return float(f"{x:.12g}")


def sturges(n):
    return math.ceil(math.log2(n)) + 1


@dataclass
class Histogram:
    edges: list          # k + 1 edges
    counts: list         # k counts
    rule: str            # how the bins were chosen
    excluded: int = 0    # values outside a user-given range
    notes: list = field(default_factory=list)

    @property
    def n(self):
        return sum(self.counts)

    @property
    def width(self):
        return self.edges[1] - self.edges[0]

    def labels(self):
        return [f"{fmt(a)}–{fmt(b)}" for a, b in zip(self.edges, self.edges[1:])]

    def percents(self):
        return [100.0 * c / self.n for c in self.counts]


def histogram(values, bins="auto", k=None, width=None, lo=None, hi=None):
    """Bin numeric values. bins: 'auto' (Sturges), 'count' (k bins) or 'width'.
    lo / hi: optional range; None uses the data's minimum / maximum."""
    if not values:
        raise HistogramError("The selection contains no numeric values.")
    user_range = lo is not None or hi is not None
    lo = min(values) if lo is None else float(lo)
    hi = max(values) if hi is None else float(hi)
    if user_range and lo >= hi:
        raise HistogramError(f"The minimum ({fmt(lo)}) must be smaller than the maximum "
                             f"({fmt(hi)}).")
    inside = [x for x in values if lo <= x <= hi]
    if not inside:
        raise HistogramError(f"No values lie between {fmt(lo)} and {fmt(hi)}.")
    excluded = len(values) - len(inside)
    n = len(inside)

    if bins == "count":
        if k is None or int(k) != k or not 1 <= k <= MAX_BINS:
            raise HistogramError(f"The number of bins must be a whole number from 1 to "
                                 f"{MAX_BINS}.")
        k = int(k)
    elif bins == "width":
        if width is None or not width > 0:
            raise HistogramError("The bin width must be a positive number.")
    elif bins != "auto":
        raise HistogramError(f"Unknown binning method {bins!r}.")

    if lo == hi:  # all values identical
        edges = [lo, lo + width] if bins == "width" else [lo, hi]
        return Histogram(edges, [n], "single bin", excluded,
                         [f"All values are equal ({fmt(lo)}), so there is a single bin."])

    if bins == "auto":
        k = sturges(n)
        rule = f"Sturges' rule: k = ceil(log2 {n}) + 1 = {k}"
    elif bins == "count":
        rule = f"{k} bins (chosen)"
    if bins == "width":
        k = max(1, math.ceil(_clean((hi - lo) / width)))
        if k > MAX_BINS:
            raise HistogramError(f"A bin width of {fmt(width)} gives {k} bins; the limit is "
                                 f"{MAX_BINS}. Choose a larger width.")
        edges = [lo] + [_clean(lo + i * width) for i in range(1, k + 1)]
        rule = f"bin width {fmt(width)} (chosen): k = ceil(({fmt(hi)} - {fmt(lo)}) / " \
               f"{fmt(width)}) = {k}"
    else:
        step = (hi - lo) / k
        edges = [lo] + [_clean(lo + i * step) for i in range(1, k)] + [hi]

    counts = [0] * k
    for x in inside:
        i = bisect.bisect_right(edges, x) - 1
        counts[min(max(i, 0), k - 1)] += 1   # the maximum goes into the last bin
    return Histogram(edges, counts, rule, excluded)
