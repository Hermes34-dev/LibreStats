"""Missing data in Calc: diagnostics report and imputation output."""

from dataclasses import dataclass

from . import calc
from . import missing as mi

IMPUTED_COLOR = 0xCCE5FF  # light blue


@dataclass
class DiagnoseOptions:
    range_text: str = ""
    has_labels: bool = True
    little: bool = True


@dataclass
class ImputeOptions:
    range_text: str = ""
    has_labels: bool = True
    method: str = "mean"          # mean | median | knn | mice
    k: int = 5
    standardize: bool = True
    mice_method: str = "pmm"      # pmm | norm
    m: int = 5
    iterations: int = 10
    seed: object = 12345          # int, or None for a random seed
    in_place: bool = False
    highlight: bool = True


@dataclass
class Table:
    """The input range: names, raw cells, and numeric variables with None for missing."""
    range: object
    names: list           # one per column of the range
    raw: list             # data rows as returned by Calc (without the label row)
    numeric: list         # indices of columns that contain at least one number
    rows: list            # numeric columns only; float or None
    first_row: int        # absolute sheet row of the first data row
    first_col: int
    text_as_missing: int  # text cells in numeric columns counted as missing

    @property
    def numeric_names(self):
        return [self.names[j] for j in self.numeric]

    @property
    def skipped_names(self):
        return [n for j, n in enumerate(self.names) if j not in self.numeric]


def read_table(doc, range_text, has_labels):
    cell_range = calc.resolve_range(doc, range_text)
    addr = cell_range.getRangeAddress()
    data = calc.used_data(cell_range)
    width = addr.EndColumn - addr.StartColumn + 1
    header = data[0] if has_labels and data else ("",) * width
    raw = [list(r) for r in data[1 if has_labels else 0:]]
    names = [str(h).strip() or f"Column {calc.column_letter(addr.StartColumn + j)}"
             for j, h in enumerate(header)]
    numeric = [j for j in range(width) if any(isinstance(r[j], float) for r in raw)]
    if not raw or not numeric:
        raise calc.RangeError("The input range contains no numbers.")
    text = sum(1 for r in raw for j in numeric if isinstance(r[j], str) and r[j].strip())
    rows = [[r[j] if isinstance(r[j], float) else None for j in numeric] for r in raw]
    return Table(cell_range, names, raw, numeric, rows,
                 addr.StartRow + (1 if has_labels else 0), addr.StartColumn, text)


def _intro(report, title, table):
    report.add(title, bold=True)
    report.add("Input range", table.range.AbsoluteName)
    report.add("Rows", len(table.rows))
    report.add("Variables analysed", ", ".join(table.numeric_names))
    if table.skipped_names:
        report.add("Skipped (no numbers)", ", ".join(table.skipped_names))
    if table.text_as_missing:
        report.add("Text cells counted as missing", table.text_as_missing)
    report.blank()


Report = calc.Report


# --- Diagnostics -------------------------------------------------------------------

def diagnose(doc, opts):
    table = read_table(doc, opts.range_text, opts.has_labels)
    rows, names = table.rows, table.numeric_names
    n, p = len(rows), len(names)
    report = Report()
    _intro(report, "LibreStats - Missing data diagnostics", table)

    report.add("Missing values by variable", bold=True)
    report.add("Variable", "n observed", "Missing", "% missing", bold=True)
    for j, name in enumerate(names):
        miss = sum(1 for r in rows if r[j] is None)
        report.add(name, n - miss, miss, 100.0 * miss / n)
    total_missing = sum(1 for r in rows for v in r if v is None)
    report.add("All variables", n * p - total_missing, total_missing,
               100.0 * total_missing / (n * p))
    complete = sum(1 for r in rows if all(v is not None for v in r))
    report.add("Complete rows", complete, n - complete, 100.0 * (n - complete) / n)
    report.blank()

    report.add("Missing-data patterns (1 = observed, 0 = missing)", bold=True)
    report.add("Rows", "% of rows", *names, "Missing in pattern", bold=True)
    for pt in mi.patterns(rows):
        report.add(pt.count, 100.0 * pt.count / n, *[int(o) for o in pt.observed],
                   sum(1 for o in pt.observed if not o))
    report.blank()

    if opts.little:
        report.add("Little's MCAR test (Little, 1988)", bold=True)
        try:
            res = mi.little_mcar(rows)
        except mi.MissingDataError as exc:
            report.add(str(exc))
        else:
            report.add("Chi-square", "df", "p", "Patterns", "Rows used", bold=True)
            report.add(res.statistic, res.df, res.p_value, res.n_patterns, res.n_rows)
            if res.p_value < 0.05:
                report.add("p < 0.05: reject MCAR. Missingness is related to the data "
                           "(MAR or MNAR); simple methods such as mean imputation or "
                           "deleting rows can bias results.")
            else:
                report.add("p >= 0.05: no evidence against MCAR. This does not prove MCAR "
                           "(the test may lack power).")
            if res.excluded_rows:
                report.add(f"{res.excluded_rows} row(s) with every value missing were excluded.")
            report.add("Note: the test assumes multivariate normality, uses EM estimates of "
                       "the mean and covariance, and cannot distinguish MAR from MNAR.")

    with calc.undo_context(doc, "LibreStats: Missing data diagnostics"):
        sheet = calc.new_sheet(doc, "Missing data")
        calc.write_rows(sheet, report.rows, report.bold)
    doc.getCurrentController().setActiveSheet(sheet)
    return sheet.getName()


# --- Imputation ----------------------------------------------------------------------

def _completed_raw(table, completed):
    """Raw rows with numeric columns replaced by the completed values."""
    out = []
    for raw_row, filled in zip(table.raw, completed):
        row = list(raw_row)
        for k, j in enumerate(table.numeric):
            row[j] = filled[k]
        out.append(row)
    return out


def impute(doc, opts):
    table = read_table(doc, opts.range_text, opts.has_labels)
    if not any(v is None for r in table.rows for v in r):
        raise calc.RangeError("There are no missing values in the numeric columns.")
    if opts.method == "mice":
        return _impute_mice(doc, opts, table)

    if opts.method in ("mean", "median"):
        res = mi.impute_simple(table.rows, opts.method)
    elif opts.method == "knn":
        res = mi.impute_knn(table.rows, opts.k, opts.standardize)
    else:
        raise ValueError(f"Unknown imputation method {opts.method!r}.")

    completed = _completed_raw(table, res.data)
    summary = Report()
    summary.add("LibreStats - " + res.method, bold=True)
    summary.add("Input range", table.range.AbsoluteName)
    if opts.method == "knn":
        summary.add("Distance", "Euclidean on " + ("z-scores" if opts.standardize else "raw values")
                    + ", over variables observed in both rows")
    summary.add("Variable", "Imputed cells", *(["Imputed value"] if res.summary else []), bold=True)
    for k, name in enumerate(table.numeric_names):
        count = sum(1 for _, j in res.imputed if j == k)
        fill = [res.summary[k][2]] if res.summary else []
        summary.add(name, count, *fill)
    for note in res.notes:
        summary.add("Note: " + note)
    if table.text_as_missing:
        summary.add(f"Note: {table.text_as_missing} text cell(s) in numeric columns were "
                    "treated as missing and replaced.")

    title = f"LibreStats: {res.method}"
    with calc.undo_context(doc, title):
        if opts.in_place:
            sheet = table.range.getSpreadsheet()
            cells = [(table.first_col + table.numeric[j], table.first_row + i)
                     for i, j in res.imputed]
            for (col, row), (i, j) in zip(cells, res.imputed):
                sheet.getCellByPosition(col, row).setValue(res.data[i][j])
            if opts.highlight:
                calc.highlight_cells(doc, sheet, cells, IMPUTED_COLOR)
            report_sheet = calc.new_sheet(doc, "Imputation")
            calc.write_rows(report_sheet, summary.rows, summary.bold)
        else:
            report_sheet = calc.new_sheet(doc, "Imputed")
            calc.write_beside(report_sheet, [table.names] + completed, summary)
            if opts.highlight:
                calc.highlight_cells(doc, report_sheet,
                                     [(table.numeric[j], i + 1) for i, j in res.imputed],
                                     IMPUTED_COLOR)
    doc.getCurrentController().setActiveSheet(report_sheet)
    return report_sheet.getName()


def _impute_mice(doc, opts, table):
    if opts.in_place:
        raise calc.RangeError("Multiple imputation creates several datasets, so it always "
                              "writes to a new sheet. Choose 'New sheet' as output.")
    res = mi.mice(table.rows, m=opts.m, iterations=opts.iterations, method=opts.mice_method,
                  seed=opts.seed)
    method = ("predictive mean matching" if res.method == "pmm"
              else "Bayesian linear regression")
    rows = [["Imputation", "Row"] + table.names]
    highlights = []
    for k, ds in enumerate(res.datasets, 1):
        for i, filled in enumerate(_completed_raw(table, ds)):
            rows.append([k, table.first_row + i + 1] + filled)
        offset = 1 + (k - 1) * len(table.raw)
        highlights.extend((2 + table.numeric[j], offset + i) for i, j in res.imputed)

    info = Report()
    info.add("LibreStats - Multiple imputation (MICE)", bold=True)
    info.add("Input range", table.range.AbsoluteName)
    info.add("Method", method)
    info.add("Imputations (m)", res.m)
    info.add("Iterations", res.iterations)
    info.add("Random seed", "random" if res.seed is None else res.seed)
    info.add("Imputed cells per dataset", len(res.imputed))
    info.add("Layout", "Datasets are stacked; 'Imputation' numbers them, "
                       "'Row' is the original sheet row.")
    for note in res.notes:
        info.add("Note: " + note)

    with calc.undo_context(doc, "LibreStats: Multiple imputation"):
        sheet = calc.new_sheet(doc, "MICE")
        calc.write_beside(sheet, rows, info)
        if opts.highlight:
            calc.highlight_cells(doc, sheet, highlights, IMPUTED_COLOR)
    doc.getCurrentController().setActiveSheet(sheet)
    return sheet.getName()
