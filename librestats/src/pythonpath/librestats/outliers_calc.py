"""Run the outlier analysis on a Calc range and write the report sheet."""

from dataclasses import dataclass

from . import calc
from . import outliers as o

HIGHLIGHT_COLOR = 0xFFCC99  # light orange


@dataclass
class OutlierOptions:
    range_text: str = ""
    has_labels: bool = True
    tukey: bool = True
    tukey_k: float = 1.5
    mad: bool = True
    mad_threshold: float = 2.5
    grubbs: bool = False
    grubbs_iterate: bool = False
    dixon: bool = False
    alpha: float = 0.05
    mahalanobis: bool = False
    mahalanobis_alpha: float = 0.001
    highlight: bool = False

    def any_method(self):
        return self.tukey or self.mad or self.grubbs or self.dixon or self.mahalanobis


@dataclass
class Variable:
    name: str
    column: int          # absolute sheet column
    rows: list           # absolute sheet row of each value
    values: list         # numeric values (non-numeric cells skipped)


def read_variables(cell_range, has_labels):
    """Read each column of the range as a variable; also return complete rows."""
    addr = cell_range.getRangeAddress()
    data = calc.used_data(cell_range)
    first = 1 if has_labels else 0
    variables = []
    for j in range(addr.EndColumn - addr.StartColumn + 1):
        col = addr.StartColumn + j
        label = data[0][j] if has_labels and data else ""
        name = str(label).strip() if str(label).strip() else f"Column {calc.column_letter(col)}"
        rows, values = [], []
        for i in range(first, len(data)):
            v = data[i][j]
            if isinstance(v, float):
                rows.append(addr.StartRow + i)
                values.append(v)
        variables.append(Variable(name, col, rows, values))

    complete_rows, complete_values = [], []
    for i in range(first, len(data)):
        row = data[i]
        if row and all(isinstance(v, float) for v in row):
            complete_rows.append(addr.StartRow + i)
            complete_values.append(list(row))
    return variables, complete_rows, complete_values


def _cell(var, flag):
    return calc.cell_name(var.column, var.rows[flag.index])


class Report:
    def __init__(self):
        self.rows = []
        self.bold = []

    def add(self, *values, bold=False):
        if bold:
            self.bold.append(len(self.rows))
        self.rows.append(list(values))

    def blank(self):
        self.rows.append([])


def _univariate_section(report, title, variables, run, score_label, highlights):
    """Summary table per variable, then the list of flagged values."""
    report.add(title, bold=True)
    results = []
    for var in variables:
        try:
            results.append((var, run(var.values)))
        except o.OutlierError as exc:
            results.append((var, exc))

    header_res = next((r for _, r in results if isinstance(r, o.Result)), None)
    labels = [label for label, _ in header_res.summary] if header_res else []
    report.add("Variable", "n", *labels, "Outliers", bold=True)
    for var, res in results:
        if isinstance(res, o.Result):
            report.add(var.name, res.n, *[v for _, v in res.summary], len(res.flags))
        else:
            report.add(var.name, len(var.values), str(res))

    flagged = [(var, f) for var, res in results if isinstance(res, o.Result) for f in res.flags]
    report.add("Flagged values", bold=True)
    if flagged:
        report.add("Variable", "Cell", "Value", "Side", score_label, bold=True)
        for var, f in flagged:
            report.add(var.name, _cell(var, f), f.value, f.side, f.score)
            highlights.add((var.column, var.rows[f.index]))
    else:
        report.add("No outliers flagged.")
    notes = {n for _, res in results if isinstance(res, o.Result) for n in res.notes}
    for note in sorted(notes):
        report.add("Note: " + note)
    report.blank()


def _grubbs_section(report, variables, opts, highlights):
    mode = "repeated until no outlier" if opts.grubbs_iterate else "single outlier"
    report.add(f"Grubbs test (two-sided, alpha = {opts.alpha:g}, {mode}; assumes normality)",
               bold=True)
    report.add("Variable", "Step", "n", "Mean", "SD", "Suspect cell", "Value", "G",
               "G critical", "p", "Outlier", bold=True)
    notes = set()
    for var in variables:
        try:
            res = o.grubbs(var.values, opts.alpha, opts.grubbs_iterate)
        except o.OutlierError as exc:
            report.add(var.name, str(exc))
            continue
        for k, s in enumerate(res.steps, 1):
            cell = calc.cell_name(var.column, var.rows[s["index"]])
            report.add(var.name, k, s["n"], s["mean"], s["sd"], cell, s["value"], s["G"],
                       s["G_crit"], s["p"], s["outlier"])
        notes.update(res.notes)
        for f in res.flags:
            highlights.add((var.column, var.rows[f.index]))
    for note in sorted(notes):
        report.add("Note: " + note)
    report.blank()


def _dixon_section(report, variables, opts, highlights):
    report.add(f"Dixon test (two-sided, alpha = {opts.alpha:g}; ratio r10/r11/r21/r22 "
               "chosen by n as in Dixon 1951; assumes normality)", bold=True)
    report.add("Variable", "n", "Ratio", "Suspect cell", "Value", "r", "r critical", "p",
               "Outlier", bold=True)
    for var in variables:
        try:
            res = o.dixon(var.values, opts.alpha)
        except o.OutlierError as exc:
            report.add(var.name, len(var.values), str(exc))
            continue
        s = res.steps[0]
        cell = calc.cell_name(var.column, var.rows[s["index"]])
        report.add(var.name, s["n"], s["ratio"], cell, s["value"], s["r"], s["r_crit"], s["p"],
                   s["outlier"])
        for f in res.flags:
            highlights.add((var.column, var.rows[f.index]))
    report.blank()


def _mahalanobis_section(report, variables, rows, values, opts, highlights):
    report.add(f"Mahalanobis distance (alpha = {opts.mahalanobis_alpha:g}; "
               "D² compared with chi-square, df = number of variables)", bold=True)
    report.add("Variables: " + ", ".join(v.name for v in variables))
    report.add(f"Complete rows used: {len(rows)} (rows with any missing or "
               "non-numeric value are excluded)")
    try:
        res = o.mahalanobis(values, opts.mahalanobis_alpha)
    except o.OutlierError as exc:
        report.add(str(exc))
        report.blank()
        return
    report.add(*[label for label, _ in res.summary], bold=True)
    report.add(*[v for _, v in res.summary])
    report.add("Flagged rows", bold=True)
    if res.flags:
        report.add("Row", "D²", "p", *[v.name for v in variables], bold=True)
        for f in res.flags:
            row = rows[f.index]
            d2, p = res.scores[f.index]
            report.add(row + 1, d2, p, *f.value)
            for v in variables:
                highlights.add((v.column, row))
    else:
        report.add("No outliers flagged.")
    for note in res.notes:
        report.add("Note: " + note)
    report.blank()


def run(doc, opts):
    """Run the selected methods; returns the name of the new report sheet."""
    cell_range = calc.resolve_range(doc, opts.range_text)
    variables, complete_rows, complete_values = read_variables(cell_range, opts.has_labels)
    if not any(v.values for v in variables):
        raise calc.RangeError("The input range contains no numbers.")

    report = Report()
    report.add("LibreStats - Outlier analysis", bold=True)
    report.add("Input range", cell_range.AbsoluteName)
    report.add("First row = variable names", opts.has_labels)
    report.add("Univariate methods test each column separately; empty and text cells are ignored.")
    report.blank()

    highlights = set()
    if opts.tukey:
        _univariate_section(
            report, f"Tukey fences (Q1 - k*IQR, Q3 + k*IQR; k = {opts.tukey_k:g}; "
                    "quartiles as Calc QUARTILE)",
            variables, lambda xs: o.tukey_fences(xs, opts.tukey_k),
            "IQRs beyond quartile", highlights)
    if opts.mad:
        _univariate_section(
            report, f"Robust z-scores (|x - median| / MAD > {opts.mad_threshold:g}; "
                    "Leys et al., 2013)",
            variables, lambda xs: o.mad_zscores(xs, opts.mad_threshold),
            "Robust z", highlights)
    if opts.grubbs:
        _grubbs_section(report, variables, opts, highlights)
    if opts.dixon:
        _dixon_section(report, variables, opts, highlights)
    if opts.mahalanobis:
        _mahalanobis_section(report, variables, complete_rows, complete_values, opts,
                             highlights)

    if opts.highlight and highlights:
        report.add(f"{len(highlights)} flagged cell(s) highlighted in orange in the input data.")

    with calc.undo_context(doc, "LibreStats: Outlier analysis"):
        sheet = calc.new_sheet(doc, "Outliers")
        calc.write_rows(sheet, report.rows, report.bold)
        if opts.highlight:
            calc.highlight_cells(doc, cell_range.getSpreadsheet(), highlights, HIGHLIGHT_COLOR)
    doc.getCurrentController().setActiveSheet(sheet)
    return sheet.getName()
