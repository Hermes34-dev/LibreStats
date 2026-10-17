"""Run the assumption checks on a Calc range and write the report sheet.

Samples come either from the columns of the input range (one sample per
column) or, if a group range is given, from one column of values split by the
labels in the group range (same rows).
"""

import math
from dataclasses import dataclass

from . import assumptions as a
from . import calc, coding

REFERENCES = {
    "sw": ["Shapiro, S. S., & Wilk, M. B. (1965). Biometrika, 52(3-4), 591-611.",
           "Royston, P. (1992). Statistics and Computing, 2(3), 117-119.",
           "Royston, P. (1995). Algorithm AS R94. Applied Statistics, 44(4), 547-551."],
    "median": ["Brown, M. B., & Forsythe, A. B. (1974). Journal of the American Statistical "
               "Association, 69(346), 364-367."],
    "mean": ["Levene, H. (1960). In I. Olkin (Ed.), Contributions to Probability and "
             "Statistics (pp. 278-292). Stanford University Press."],
}


@dataclass
class NormalityOptions:
    range_text: str = ""
    has_labels: bool = True
    group_range: str = ""
    alpha: float = 0.05


@dataclass
class VarianceOptions:
    range_text: str = ""
    has_labels: bool = True
    group_range: str = ""
    center: str = "median"      # median (Brown-Forsythe) | mean (Levene)
    alpha: float = 0.05


@dataclass
class Sample:
    name: str
    values: list
    text: int = 0           # non-empty text cells skipped


@dataclass
class Samples:
    items: list
    input_name: str
    group_name: str = ""    # address of the group range, if used
    variable: str = ""      # name of the value column when grouped
    ungrouped: int = 0      # values without a group label


def _numbers(cells):
    values = [v for v in cells if isinstance(v, float)]
    text = sum(1 for v in cells if isinstance(v, str) and v.strip())
    return values, text


def read_samples(doc, range_text, has_labels, group_range=""):
    cell_range = calc.resolve_range(doc, range_text)
    addr = cell_range.getRangeAddress()
    first = 1 if has_labels else 0

    if not group_range.strip():
        data = calc.used_data(cell_range)
        if len(data) <= first:
            raise calc.RangeError("The input range contains no data rows.")
        items = []
        for j in range(addr.EndColumn - addr.StartColumn + 1):
            label = str(data[0][j]).strip() if has_labels else ""
            name = label or f"Column {calc.column_letter(addr.StartColumn + j)}"
            values, text = _numbers([r[j] for r in data[first:]])
            items.append(Sample(name, values, text))
        return Samples(items, cell_range.AbsoluteName)

    groups = calc.resolve_range(doc, group_range)
    gaddr = groups.getRangeAddress()
    if addr.StartColumn != addr.EndColumn or gaddr.StartColumn != gaddr.EndColumn:
        raise calc.RangeError("With group labels, the input range and the group range must "
                              "each be a single column.")
    values, labels = cell_range.getDataArray(), groups.getDataArray()
    if len(values) != len(labels):
        raise calc.RangeError(f"The input range has {len(values)} rows but the group range has "
                              f"{len(labels)}; they must cover the same rows.")
    rows = max(len(calc.used_data(cell_range)), len(calc.used_data(groups)))
    if rows <= first:
        raise calc.RangeError("The input range contains no data rows.")
    variable = str(values[0][0]).strip() if has_labels else ""
    by_group, ungrouped, text = {}, 0, {}
    for (v,), (g,) in zip(values[first:rows], labels[first:rows]):
        if coding.is_missing(g):
            ungrouped += isinstance(v, float)
            continue
        key = coding.label(g)
        by_group.setdefault(key, [])
        if isinstance(v, float):
            by_group[key].append(v)
        elif isinstance(v, str) and v.strip():
            text[key] = text.get(key, 0) + 1
    if not by_group:
        raise calc.RangeError("The group range contains no group labels.")
    items = [Sample(k, vals, text.get(k, 0)) for k, vals in by_group.items()]
    return Samples(items, cell_range.AbsoluteName, groups.AbsoluteName, variable, ungrouped)


def _mean_sd(xs):
    n = len(xs)
    mean = sum(xs) / n if n else None
    sd = math.sqrt(sum((x - mean) ** 2 for x in xs) / (n - 1)) if n > 1 else None
    return mean, sd


def _conclusion(p, alpha, rejected, kept):
    return f"{rejected} (p < {alpha:g})" if p < alpha else kept


def _intro(report, title, samples, settings):
    report.add("LibreStats - " + title, bold=True)
    report.add("Input range", samples.input_name)
    if samples.group_name:
        report.add("Group labels", samples.group_name)
        if samples.variable:
            report.add("Variable", samples.variable)
    for label, value in settings:
        report.add(label, value)
    report.blank()


def _notes(report, samples, notes, refs):
    if samples.ungrouped:
        notes.append(f"{samples.ungrouped} value(s) without a group label were not used.")
    text = sum(s.text for s in samples.items)
    if text:
        notes.append(f"{text} text cell(s) ignored.")
    if notes:
        report.blank()
        for n in dict.fromkeys(notes):
            report.add("Note: " + n)
    report.blank()
    report.add("References", bold=True)
    for r in refs:
        report.add(r)


def _write(doc, report, sheet_name, undo_title):
    with calc.undo_context(doc, "LibreStats: " + undo_title):
        sheet = calc.new_sheet(doc, sheet_name)
        calc.write_rows(sheet, report.rows, report.bold)
    doc.getCurrentController().setActiveSheet(sheet)
    return sheet.getName()


def run_normality(doc, opts):
    samples = read_samples(doc, opts.range_text, opts.has_labels, opts.group_range)
    report = calc.Report()
    _intro(report, "Normality test (Shapiro-Wilk)", samples,
           [("Null hypothesis", "the data come from a normal distribution"),
            ("Significance level", opts.alpha)])
    report.add("Group" if samples.group_name else "Variable", "n", "Mean", "SD", "W",
               "p-value", "Conclusion", bold=True)
    notes, computed = [], 0
    for s in samples.items:
        mean, sd = _mean_sd(s.values)
        try:
            res = a.shapiro_wilk(s.values)
        except a.AssumptionError as exc:
            report.add(s.name, len(s.values), mean, sd, "not computed: " + str(exc))
            continue
        computed += 1
        report.add(s.name, res.n, mean, sd, res.statistic, res.p,
                   _conclusion(res.p, opts.alpha, "normality rejected",
                               "no evidence against normality"))
        notes += res.notes
    if not computed:
        raise calc.RangeError("Shapiro-Wilk could not be computed for any variable "
                              "(each needs at least 3 values that are not all equal).")
    notes.append("With large samples the test detects deviations from normality too small to "
                 "matter; with small samples it may miss real ones. Judge it together with a "
                 "histogram or Q-Q plot.")
    _notes(report, samples, notes, REFERENCES["sw"])
    return _write(doc, report, "Normality", "Shapiro-Wilk test")


def run_variances(doc, opts):
    samples = read_samples(doc, opts.range_text, opts.has_labels, opts.group_range)
    title = ("Equal variances: Levene's test, Brown-Forsythe version (median-centred)"
             if opts.center == "median" else "Equal variances: Levene's test (mean-centred)")
    res = a.levene([s.values for s in samples.items], opts.center)
    report = calc.Report()
    _intro(report, title, samples,
           [("Null hypothesis", "all groups have the same variance"),
            ("Centre", "group median" if opts.center == "median" else "group mean"),
            ("Significance level", opts.alpha)])
    report.add("Group", "n", "Mean", "Median", "SD", "Variance", bold=True)
    variances = []
    for s in samples.items:
        mean, sd = _mean_sd(s.values)
        variances.append(sd * sd)
        report.add(s.name, len(s.values), mean, a._median(s.values), sd, sd * sd)
    report.blank()
    df = dict(res.params)
    report.add("F", "df1", "df2", "p-value", "Conclusion", bold=True)
    report.add(res.statistic, df["df1"], df["df2"], res.p,
               _conclusion(res.p, opts.alpha, "variances differ",
                           "no evidence of unequal variances"))
    report.blank()
    smallest = min(variances)
    report.add("Largest / smallest variance", max(variances) / smallest if smallest else "")
    notes = ["The test uses the absolute deviations of each value from its group "
             f"{opts.center}, compared by one-way ANOVA."]
    if opts.center == "median":
        notes.append("The median-centred version keeps its error rate when the data are "
                     "skewed or heavy-tailed (Brown & Forsythe, 1974); the mean-centred "
                     "original is only accurate for normal data.")
    _notes(report, samples, notes, REFERENCES[opts.center])
    return _write(doc, report, "Equal variances", "Equal variances test")


RUNNERS = {"normality": run_normality, "variances": run_variances}
