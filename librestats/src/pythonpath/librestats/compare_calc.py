"""Run the comparison tests on a Calc range and write the report sheet.

Independent samples (one-sample and Welch t-test, one-way and Welch ANOVA) are
read as in the assumption checks: one sample per column, or one column of
values split by a group range. Paired samples (paired t-test, repeated-measures
ANOVA) are columns of the input range, one row per subject; rows with a
missing value are left out.
"""

import math
from dataclasses import dataclass

from . import calc
from . import compare as c
from .assumptions_calc import (Samples, _conclusion, _intro, _mean_sd, _notes, _numbers, _write,
                               read_samples)

REFERENCES = {
    "t": ["Student (1908). The probable error of a mean. Biometrika, 6(1), 1-25."],
    "welch_t": ["Welch, B. L. (1947). Biometrika, 34(1-2), 28-35.",
                "Satterthwaite, F. E. (1946). Biometrics Bulletin, 2(6), 110-114.",
                "Delacre, M., Lakens, D., & Leys, C. (2017). International Review of Social "
                "Psychology, 30(1), 92-101.",
                "Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences "
                "(2nd ed.). Erlbaum."],
    "d": ["Lakens, D. (2013). Frontiers in Psychology, 4, 863."],
    "anova": ["Fisher, R. A. (1925). Statistical Methods for Research Workers. Oliver & Boyd.",
              "Olejnik, S., & Algina, J. (2003). Psychological Methods, 8(4), 434-447."],
    "welch_anova": ["Welch, B. L. (1951). Biometrika, 38(3-4), 330-336.",
                    "Delacre, M., Leys, C., Mora, Y. L., & Lakens, D. (2019). International "
                    "Review of Social Psychology, 32(1), 13."],
    "rm": ["Mauchly, J. W. (1940). Annals of Mathematical Statistics, 11(2), 204-209.",
           "Greenhouse, S. W., & Geisser, S. (1959). Psychometrika, 24(2), 95-112.",
           "Huynh, H., & Feldt, L. S. (1976). Journal of Educational Statistics, 1(1), 69-82.",
           "Lecoutre, B. (1991). Journal of Educational Statistics, 16(4), 371-372.",
           "Girden, E. R. (1992). ANOVA: Repeated Measures. Sage.",
           "Bakeman, R. (2005). Behavior Research Methods, 37(3), 379-384."],
}

ALTERNATIVE_SIGNS = {"two-sided": "≠", "greater": ">", "less": "<"}


@dataclass
class TTestOptions:
    range_text: str = ""
    has_labels: bool = True
    group_range: str = ""
    mu: float = 0.0                 # one-sample test only
    alternative: str = "two-sided"  # two-sided | greater | less
    alpha: float = 0.05


@dataclass
class AnovaOptions:
    range_text: str = ""
    has_labels: bool = True
    group_range: str = ""
    alpha: float = 0.05


@dataclass
class Paired:
    names: list
    rows: list              # complete rows only
    input_name: str
    incomplete: int = 0     # rows with a missing value in some column
    text: int = 0


def read_paired(doc, range_text, has_labels, minimum=2, maximum=None):
    """Columns of the input range as repeated measurements, one row per subject."""
    cell_range = calc.resolve_range(doc, range_text)
    addr = cell_range.getRangeAddress()
    k = addr.EndColumn - addr.StartColumn + 1
    if k < minimum or (maximum and k > maximum):
        need = f"exactly {minimum}" if minimum == maximum else f"at least {minimum}"
        raise calc.RangeError(f"Select {need} columns, one per condition ({k} selected).")
    data = calc.used_data(cell_range)
    first = 1 if has_labels else 0
    if len(data) <= first:
        raise calc.RangeError("The input range contains no data rows.")
    names = [(str(data[0][j]).strip() if has_labels else "") or
             f"Column {calc.column_letter(addr.StartColumn + j)}" for j in range(k)]
    rows, incomplete, text = [], 0, 0
    for r in data[first:]:
        text += _numbers(r)[1]
        if all(isinstance(v, float) for v in r):
            rows.append(list(r))
        elif any(isinstance(v, float) for v in r):
            incomplete += 1
    return Paired(names, rows, cell_range.AbsoluteName, incomplete, text)


def _bound(v):
    return "∞" if v == math.inf else "-∞" if v == -math.inf else v


def _confidence(alpha):
    return f"{100 * (1 - alpha):g}%"


def _descriptives(report, heading, items):
    """items: (name, values). Rows: name, n, mean, SD, SE."""
    report.add(heading, "n", "Mean", "SD", "SE", bold=True)
    for name, values in items:
        mean, sd = _mean_sd(values)
        report.add(name, len(values), mean, sd, sd / math.sqrt(len(values)) if sd else "")
    report.blank()


def _paired_notes(report, paired, notes, refs):
    if paired.incomplete:
        notes.append(f"{paired.incomplete} row(s) with a missing value in some column were "
                     "left out (only complete rows are used).")
    if paired.text:
        notes.append(f"{paired.text} text cell(s) ignored.")
    _notes(report, Samples([], paired.input_name), notes, refs)


def _effects(report, effects):
    report.blank()
    report.add("Effect size", "Value", bold=True)
    for label, value in effects:
        report.add(label, value)


def _anova_table(report, res, alpha):
    report.add("Source", "SS", "df", "MS", "F", "p-value", bold=True)
    for source, ss, df, ms, f, p in res.table:
        report.add(source, ss, df, "" if ms is None else ms, "" if f is None else f,
                   "" if p is None else p)
    report.blank()
    report.add("Conclusion", _conclusion(res.p, alpha, "means differ",
                                          "no evidence that the means differ"))


# --- t-tests ----------------------------------------------------------------------------

def run_one_sample(doc, opts):
    samples = read_samples(doc, opts.range_text, opts.has_labels, opts.group_range)
    sign = ALTERNATIVE_SIGNS[opts.alternative]
    report = calc.Report()
    _intro(report, "One-sample t-test", samples,
           [("Test value (μ0)", opts.mu),
            ("Alternative hypothesis", f"mean {sign} {opts.mu:g}"),
            ("Confidence level", _confidence(opts.alpha)),
            ("Significance level", opts.alpha)])
    conf = _confidence(opts.alpha)
    report.add("Group" if samples.group_name else "Variable", "n", "Mean", "SD", "SE",
               "Mean - μ0", "t", "df", "p-value", f"{conf} CI lower", f"{conf} CI upper",
               "Cohen's d", "Conclusion", bold=True)
    computed = 0
    for s in samples.items:
        try:
            res = c.one_sample_t(s.values, opts.mu, opts.alternative, opts.alpha)
        except c.CompareError as exc:
            report.add(s.name, len(s.values), "not computed: " + str(exc))
            continue
        computed += 1
        report.add(s.name, res.n, res.estimate, res.se * math.sqrt(res.n), res.se,
                   res.estimate - opts.mu, res.statistic, res.df, res.p, _bound(res.ci[0]),
                   _bound(res.ci[1]), res.effect,
                   _conclusion(res.p, opts.alpha, f"mean {sign} {opts.mu:g}",
                                "no evidence against H0"))
    if not computed:
        raise calc.RangeError("The t-test could not be computed for any variable (each needs "
                              "at least 2 values that are not all equal).")
    notes = ["The confidence interval is for the mean.",
             "Cohen's d = (mean - μ0) / SD.",
             "The test assumes independent values from a normal distribution; with larger "
             "samples it is robust to moderate non-normality."]
    _notes(report, samples, notes, REFERENCES["t"] + REFERENCES["d"])
    return _write(doc, report, "One-sample t-test", "One-sample t-test")


def _two_samples(samples):
    if len(samples.items) != 2:
        where = "groups in the group range" if samples.group_name else "columns"
        raise calc.RangeError(f"Welch's t-test compares exactly 2 samples; the input has "
                              f"{len(samples.items)} {where}.")
    return samples.items


def run_welch(doc, opts):
    samples = read_samples(doc, opts.range_text, opts.has_labels, opts.group_range)
    s1, s2 = _two_samples(samples)
    res = c.welch_t(s1.values, s2.values, opts.alternative, opts.alpha)
    sign = ALTERNATIVE_SIGNS[opts.alternative]
    hypothesis = f"mean({s1.name}) {sign} mean({s2.name})"
    report = calc.Report()
    _intro(report, "Welch's t-test (unequal variances)", samples,
           [("Alternative hypothesis", hypothesis),
            ("Confidence level", _confidence(opts.alpha)),
            ("Significance level", opts.alpha)])
    _descriptives(report, "Group", [(s1.name, s1.values), (s2.name, s2.values)])
    conf = _confidence(opts.alpha)
    report.add("Difference", "SE", "t", "df", "p-value", f"{conf} CI lower", f"{conf} CI upper",
               "Cohen's d", "Conclusion", bold=True)
    report.add(res.estimate, res.se, res.statistic, res.df, res.p, _bound(res.ci[0]),
               _bound(res.ci[1]), res.effect,
               _conclusion(res.p, opts.alpha, hypothesis, "no evidence that the means differ"))
    notes = [f"Difference = mean({s1.name}) - mean({s2.name}); the confidence interval is for "
             "this difference.",
             "Welch's test does not assume equal variances, and is recommended over Student's "
             "pooled t-test as the default (Delacre et al., 2017). The degrees of freedom are "
             "the Welch-Satterthwaite approximation.",
             "Cohen's d uses the square root of the average of the two variances as "
             "standardiser (Cohen, 1988), so that it does not assume equal variances either."]
    _notes(report, samples, notes, REFERENCES["welch_t"])
    return _write(doc, report, "Welch t-test", "Welch's t-test")


def _correlation(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else ""


def run_paired(doc, opts):
    paired = read_paired(doc, opts.range_text, opts.has_labels, 2, 2)
    a, b = paired.names
    x1, x2 = [r[0] for r in paired.rows], [r[1] for r in paired.rows]
    res = c.paired_t(x1, x2, opts.alternative, opts.alpha)
    sign = ALTERNATIVE_SIGNS[opts.alternative]
    report = calc.Report()
    _intro(report, "Paired t-test", Samples([], paired.input_name),
           [("Alternative hypothesis", f"mean difference ({a} - {b}) {sign} 0"),
            ("Confidence level", _confidence(opts.alpha)),
            ("Significance level", opts.alpha)])
    _descriptives(report, "Variable", [(a, x1), (b, x2),
                                       (f"Difference ({a} - {b})",
                                        [u - v for u, v in zip(x1, x2)])])
    report.add("Correlation of pairs (r)", _correlation(x1, x2))
    report.blank()
    conf = _confidence(opts.alpha)
    report.add("Mean difference", "SE", "t", "df", "p-value", f"{conf} CI lower",
               f"{conf} CI upper", "Cohen's d_z", "Conclusion", bold=True)
    report.add(res.estimate, res.se, res.statistic, res.df, res.p, _bound(res.ci[0]),
               _bound(res.ci[1]), res.effect,
               _conclusion(res.p, opts.alpha, f"mean difference {sign} 0",
                            "no evidence of a difference"))
    notes = ["Each row is one pair. The test is a one-sample t-test of the differences; it "
             "assumes they are independent and normally distributed.",
             "Cohen's d_z = mean difference / SD of the differences (Lakens, 2013)."]
    _paired_notes(report, paired, notes, REFERENCES["t"] + REFERENCES["d"])
    return _write(doc, report, "Paired t-test", "Paired t-test")


# --- ANOVA ------------------------------------------------------------------------------

def _anova_samples(samples):
    if len(samples.items) < 2:
        raise calc.RangeError("The ANOVA compares at least 2 groups: select several columns, "
                              "or one column of values with a group range.")
    return samples.items


def run_anova(doc, opts):
    samples = read_samples(doc, opts.range_text, opts.has_labels, opts.group_range)
    items = _anova_samples(samples)
    res = c.one_way_anova([s.values for s in items])
    report = calc.Report()
    _intro(report, "One-way ANOVA", samples,
           [("Null hypothesis", "all group means are equal"),
            ("Significance level", opts.alpha)])
    _descriptives(report, "Group", [(s.name, s.values) for s in items])
    _anova_table(report, res, opts.alpha)
    _effects(report, res.effects)
    notes = ["The classical F-test assumes normal data and equal variances in all groups. "
             "With unequal variances, especially with unequal group sizes, its error rate is "
             "wrong; Welch's ANOVA does not need equal variances (Delacre et al., 2019).",
             "ω² is less biased than η² as an estimate of the population effect "
             "(Olejnik & Algina, 2003)."]
    _notes(report, samples, notes, REFERENCES["anova"])
    return _write(doc, report, "One-way ANOVA", "One-way ANOVA")


def run_welch_anova(doc, opts):
    samples = read_samples(doc, opts.range_text, opts.has_labels, opts.group_range)
    items = _anova_samples(samples)
    res = c.welch_anova([s.values for s in items])
    report = calc.Report()
    _intro(report, "Welch's ANOVA (unequal variances)", samples,
           [("Null hypothesis", "all group means are equal"),
            ("Significance level", opts.alpha)])
    _descriptives(report, "Group", [(s.name, s.values) for s in items])
    report.add("F", "df1", "df2", "p-value", "Conclusion", bold=True)
    report.add(res.statistic, res.df1, res.df2, res.p,
               _conclusion(res.p, opts.alpha, "means differ",
                            "no evidence that the means differ"))
    notes = ["Welch's F weights each group by n / variance and does not assume equal "
             "variances; it assumes normal data within each group. df2 is not an integer.",
             "With 2 groups it equals Welch's t-test (F = t²)."]
    _notes(report, samples, notes, REFERENCES["welch_anova"])
    return _write(doc, report, "Welch ANOVA", "Welch's ANOVA")


def run_rm_anova(doc, opts):
    paired = read_paired(doc, opts.range_text, opts.has_labels, 2)
    res = c.rm_anova(paired.rows)
    k = len(paired.names)
    report = calc.Report()
    _intro(report, "Repeated-measures ANOVA (one within-subjects factor)",
           Samples([], paired.input_name),
           [("Null hypothesis", "all condition means are equal"),
            ("Number of subjects (complete rows)", res.n),
            ("Number of conditions", k),
            ("Significance level", opts.alpha)])
    _descriptives(report, "Condition",
                  [(name, [r[j] for r in paired.rows]) for j, name in enumerate(paired.names)])
    _anova_table(report, res, opts.alpha)
    sph = res.sphericity
    if "W" in sph:
        report.blank()
        report.add("Mauchly's test of sphericity", bold=True)
        report.add("W", "χ²", "df", "p-value", "Conclusion", bold=True)
        report.add(sph["W"], sph["chi2"], sph["df"], sph["p"],
                   _conclusion(sph["p"], opts.alpha, "sphericity rejected",
                                "no evidence against sphericity"))
    report.blank()
    report.add("Correction", "Epsilon (ε)", "df1", "df2", "p-value", "Conclusion", bold=True)
    for name, eps, df1, df2, p in res.corrections:
        report.add(name, eps, df1, df2, p,
                   _conclusion(p, opts.alpha, "means differ",
                                "no evidence that the means differ"))
    _effects(report, res.effects)
    notes = list(res.notes)
    if k > 2:
        notes.append("Mauchly's test has little power in small samples and too much in large "
                     "ones, so many authors report a corrected p-value regardless. A common "
                     "rule (Girden, 1992): use Huynh-Feldt when the Greenhouse-Geisser ε > 0.75, "
                     "otherwise Greenhouse-Geisser.")
    notes.append("Generalized η² is comparable across within- and between-subjects designs "
                 "(Olejnik & Algina, 2003; Bakeman, 2005).")
    _paired_notes(report, paired, notes, REFERENCES["rm"] + REFERENCES["anova"][1:])
    return _write(doc, report, "RM ANOVA", "Repeated-measures ANOVA")


RUNNERS = {"ttestOne": run_one_sample, "ttestWelch": run_welch, "ttestPaired": run_paired,
           "anovaOneWay": run_anova, "anovaWelch": run_welch_anova, "anovaRM": run_rm_anova}
