"""Normal quantile-quantile (Q-Q) plot coordinates and reference lines.

The i-th smallest value is plotted against the standard normal quantile
Phi^-1((i - a) / (n + 1 - 2a)), with a the plotting-position constant
(Blom a = 3/8, van der Waerden 0, Tukey 1/3, Rankit 1/2). Tied values get
consecutive positions, as in R's qqnorm.
"""

from dataclasses import dataclass, field

from .distributions import norm_ppf
from .outliers import quantile
from .transform import NORMAL_SCORE_NAMES, NORMAL_SCORES

LINES = ("quartiles", "meansd", "none")


class QQError(ValueError):
    pass


@dataclass
class QQ:
    theoretical: list
    sample: list                 # sorted values
    intercept: float = None      # reference line y = intercept + slope * z (None: no line)
    slope: float = None
    notes: list = field(default_factory=list)

    def line(self):
        if self.slope is None:
            return []
        return [self.intercept + self.slope * z for z in self.theoretical]


def plotting_positions(n, method="blom"):
    a = NORMAL_SCORES[method]
    return [norm_ppf((i - a) / (n + 1 - 2 * a)) for i in range(1, n + 1)]


def qq_normal(values, method="blom", line="quartiles"):
    if method not in NORMAL_SCORES:
        raise QQError(f"Unknown plotting positions {method!r}.")
    if line not in LINES:
        raise QQError(f"Unknown reference line {line!r}.")
    x = sorted(v for v in values if v is not None)
    n = len(x)
    if n < 3:
        raise QQError("A Q-Q plot needs at least 3 values.")
    res = QQ(plotting_positions(n, method), x)
    if line == "quartiles":
        # through the first and third quartiles, as R's qqline (quantile type 7)
        q1, q3 = quantile(x, 0.25), quantile(x, 0.75)
        z3 = norm_ppf(0.75)
        res.slope = (q3 - q1) / (2 * z3)
        res.intercept = (q1 + q3) / 2
    elif line == "meansd":
        mean = sum(x) / n
        res.slope = (sum((v - mean) ** 2 for v in x) / (n - 1)) ** 0.5
        res.intercept = mean
    if x[0] == x[-1]:
        res.notes.append("All values are equal.")
    return res


def describe(method):
    a = NORMAL_SCORES[method]
    return f"{NORMAL_SCORE_NAMES[method]}: (i - {a:.4g}) / (n + {1 - 2 * a:.4g})"
