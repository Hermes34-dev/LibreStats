"""Transformations, scaling, log-ratios and (re)coding on a Calc range.

Every tool adds new columns, either on a new sheet (copy of the input plus the
new columns, with a report beside it) or directly to the right of the input
range (report on its own sheet). Settings come in PrepOptions.settings:

  transform  method=log|sqrt|boxcox|yeojohnson|logit, base=e|10|2,
             offset=none|one|half_min|<number>, lambda=<number or "">,
             percent=bool, squeeze=bool
  scale      method=centre|z|range|pareto
  logratio   method=clr|ilr, zeros=error|replace
  dummies    drop_first=bool
  presence   threshold=<number>, empty_absent=bool
  bins       method=width|quantile|kmeans, k=int
  recode     rules=<text, see coding.parse_rules>
  ranks      ranks=bool, ties=average|min|max|dense, scores=bool,
             score_method=blom|vdw|tukey|rankit
"""

from dataclasses import dataclass, field

from . import calc, coding, compositional, transform
from .coding import CodingError
from .compositional import CompositionError
from .transform import TransformError

METHOD_ERRORS = (TransformError, CodingError, CompositionError)


@dataclass
class PrepOptions:
    range_text: str = ""
    has_labels: bool = True
    beside: bool = False
    settings: dict = field(default_factory=dict)


@dataclass
class Input:
    range: object
    names: list
    raw: list            # data rows (without the label row)
    has_labels: bool
    first_row: int       # sheet row of the label row (or first data row)
    first_col: int
    last_col: int

    def numeric(self, j):
        """Column j as floats with None for empty/text cells; also the text count."""
        vals, text = [], 0
        for r in self.raw:
            v = r[j]
            if isinstance(v, float):
                vals.append(v)
            else:
                vals.append(None)
                text += isinstance(v, str) and bool(v.strip())
        return vals, text

    def column(self, j):
        return [None if coding.is_missing(r[j]) else r[j] for r in self.raw]


def read_input(doc, opts):
    cell_range = calc.resolve_range(doc, opts.range_text)
    a = cell_range.getRangeAddress()
    data = calc.used_data(cell_range)
    if len(data) < (2 if opts.has_labels else 1):
        raise calc.RangeError("The input range contains no data rows.")
    width = a.EndColumn - a.StartColumn + 1
    header = data[0] if opts.has_labels else ("",) * width
    names = [calc.column_letter(a.StartColumn + j) if not str(h).strip() else coding.label(h)
             for j, h in enumerate(header)]
    raw = [list(r) for r in data[1 if opts.has_labels else 0:]]
    return Input(cell_range, names, raw, opts.has_labels, a.StartRow, a.StartColumn, a.EndColumn)


# --- Generic per-variable runner -----------------------------------------------------

@dataclass
class Outcome:
    """Result for one input variable: new columns plus report material."""
    variable: str
    columns: list = field(default_factory=list)   # [(header, values)]
    params: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    table: list = field(default_factory=list)
    error: str = ""


def _run_each(inp, func, numeric=True):
    """Apply func(values) -> Outcome-like tuple to every column."""
    outcomes = []
    for j, name in enumerate(inp.names):
        if numeric:
            values, text = inp.numeric(j)
            if not any(v is not None for v in values):
                outcomes.append(Outcome(name, error="no numeric values"))
                continue
        else:
            values, text = inp.column(j), 0
        try:
            columns, params, notes, table = func(name, values)
        except METHOD_ERRORS as exc:
            outcomes.append(Outcome(name, error=str(exc)))
            continue
        if text:
            notes = list(notes) + [f"{text} text cell(s) treated as missing."]
        outcomes.append(Outcome(name, columns, params, notes, table))
    if not any(o.columns for o in outcomes):
        problems = "; ".join(f"{o.variable}: {o.error}" for o in outcomes)
        raise calc.RangeError("Nothing could be computed. " + problems)
    return outcomes


def _report(title, inp, settings, outcomes):
    rep = calc.Report()
    rep.add("LibreStats - " + title, bold=True)
    rep.add("Input range", inp.range.AbsoluteName)
    for label, value in settings:
        rep.add(label, value)
    rep.blank()
    labels = next((o.params for o in outcomes if o.columns), [])
    rep.add("Variable", "New column(s)", *[lab for lab, _ in labels], bold=True)
    for o in outcomes:
        if o.error:
            rep.add(o.variable, "not computed: " + o.error)
        else:
            rep.add(o.variable, ", ".join(h for h, _ in o.columns), *[v for _, v in o.params])
    for o in outcomes:
        if o.table:
            rep.blank()
            rep.add(o.variable, bold=True)
            for i, row in enumerate(o.table):
                rep.add(*row, bold=i == 0)
    notes = []
    for o in outcomes:
        for n in o.notes:
            if n not in notes:
                notes.append(n)
    if notes:
        rep.blank()
        for n in notes:
            rep.add("Note: " + n)
    return rep


def _write(doc, inp, outcomes, report, opts, sheet_name, undo_title):
    new = [c for o in outcomes for c in o.columns]
    headers = [h for h, _ in new]
    n = len(inp.raw)
    values = [[col[i] for _, col in new] for i in range(n)]
    with calc.undo_context(doc, "LibreStats: " + undo_title):
        if opts.beside:
            sheet = inp.range.getSpreadsheet()
            top = inp.first_row
            rows = ([headers] if inp.has_labels else []) + values
            target = sheet.getCellRangeByPosition(inp.last_col + 1, top,
                                                  inp.last_col + len(headers), top + len(rows) - 1)
            if target.queryContentCells(1 | 2 | 4 | 16).getCount():
                raise calc.RangeError(f"The {len(headers)} column(s) to the right of the input "
                                      "range are not empty. Clear them or choose 'New sheet'.")
            target.setDataArray(tuple(tuple(calc._cell_value(v) for v in r) for r in rows))
            if inp.has_labels:
                sheet.getCellRangeByPosition(inp.last_col + 1, top, inp.last_col + len(headers),
                                             top).CharWeight = calc.BOLD
            report_sheet = calc.new_sheet(doc, sheet_name + " report")
            calc.write_rows(report_sheet, report.rows, report.bold)
        else:
            report_sheet = calc.new_sheet(doc, sheet_name)
            table = [inp.names + headers] + [r + v for r, v in zip(inp.raw, values)]
            calc.write_beside(report_sheet, table, report)
    doc.getCurrentController().setActiveSheet(report_sheet)
    return report_sheet.getName()


def _flag(settings, key, default=False):
    v = settings.get(key, default)
    return v if isinstance(v, bool) else str(v).lower() in ("1", "true", "yes")


# --- Transformations ----------------------------------------------------------------

def run_transform(doc, opts):
    s = opts.settings
    method = s.get("method", "log")
    base = str(s.get("base", "e"))
    offset = s.get("offset", "none")
    if offset not in transform.OFFSETS:
        offset = float(offset)
    lam = s.get("lambda", "")
    lam = None if lam in ("", None) else float(lam)
    percent, squeeze = _flag(s, "percent"), _flag(s, "squeeze", True)

    def off_text(res):
        c = dict(res.params).get("Offset added")
        return f"+{c:g}" if c else ""

    def func(name, values):
        if method == "log":
            res = transform.log_transform(values, base, offset)
            fn = {"e": "ln", "10": "log10", "2": "log2"}[base]
            header = f"{fn}({name}{off_text(res)})"
        elif method == "sqrt":
            res = transform.sqrt_transform(values, offset)
            header = f"sqrt({name}{off_text(res)})"
        elif method == "boxcox":
            res = transform.boxcox(values, lam, offset)
            header = f"BoxCox({name}{off_text(res)})"
        elif method == "yeojohnson":
            res = transform.yeojohnson(values, lam)
            header = f"YeoJohnson({name})"
        elif method == "logit":
            res = transform.logit_transform(values, percent, squeeze)
            header = f"logit({name})"
        else:
            raise ValueError(f"Unknown transformation {method!r}")
        return [(header, res.values)], res.params, res.notes, []

    titles = {"log": "Log transformation", "sqrt": "Square root transformation",
              "boxcox": "Box-Cox transformation (Box & Cox, 1964)",
              "yeojohnson": "Yeo-Johnson transformation (Yeo & Johnson, 2000)",
              "logit": "Logit transformation"}
    inp = read_input(doc, opts)
    outcomes = _run_each(inp, func)
    settings = []
    if method in ("boxcox", "yeojohnson"):
        settings.append(("Lambda", "maximum likelihood, 95% likelihood-ratio CI"
                         if lam is None else lam))
    report = _report(titles[method], inp, settings, outcomes)
    return _write(doc, inp, outcomes, report, opts, "Transformed", titles[method])


# --- Scaling ------------------------------------------------------------------------

SCALE_TITLES = {"centre": ("Centring", "centred"), "z": ("Autoscaling (z-scores)", "z"),
                "range": ("Range scaling (0-1)", "range01"), "pareto": ("Pareto scaling", "pareto")}


def run_scale(doc, opts):
    method = opts.settings.get("method", "z")
    title, prefix = SCALE_TITLES[method]

    def func(name, values):
        res = transform.scale(values, method)
        return [(f"{prefix}({name})", res.values)], res.params, res.notes, []

    inp = read_input(doc, opts)
    outcomes = _run_each(inp, func)
    report = _report(title, inp, [("Reference", "van den Berg et al. (2006)")], outcomes)
    return _write(doc, inp, outcomes, report, opts, "Scaled", title)


# --- Compositional data -------------------------------------------------------------

def run_logratio(doc, opts):
    s = opts.settings
    method = s.get("method", "clr")
    inp = read_input(doc, opts)
    cols, text = [], 0
    for j in range(len(inp.names)):
        vals, t = inp.numeric(j)
        if not any(v is not None for v in vals):
            raise calc.RangeError(f"Column '{inp.names[j]}' has no numbers. Select only the "
                                  "columns that are parts of the composition.")
        cols.append(vals)
        text += t
    rows = [[cols[j][i] for j in range(len(cols))] for i in range(len(inp.raw))]
    res = compositional.log_ratios(rows, inp.names, method,
                                   "replace" if s.get("zeros") == "replace" else "error")
    new = [(name, [None if r is None else r[k] for r in res.rows])
           for k, name in enumerate(res.names)]
    notes = list(res.notes) + ([f"{text} text cell(s) treated as missing."] if text else [])
    outcome = Outcome("All parts", new, res.params, notes)
    title = "Log-ratio transformation (Aitchison, 1982)"
    report = _report(title, inp, [], [outcome])
    return _write(doc, inp, [outcome], report, opts, "Log-ratios", "Log-ratio transformation")


# --- Coding -------------------------------------------------------------------------

def _coded(name, res, fmt):
    return ([(fmt.format(name=name, suffix=s), col) for s, col in zip(res.names, res.columns)],
            res.params, res.notes, res.table)


def run_dummies(doc, opts):
    drop = _flag(opts.settings, "drop_first")
    inp = read_input(doc, opts)
    outcomes = _run_each(inp, lambda name, v: _coded(name, coding.dummies(v, drop),
                                                      "{name}={suffix}"), numeric=False)
    report = _report("Dummy variables (one-hot coding)", inp,
                     [("Reference category dropped", drop)], outcomes)
    return _write(doc, inp, outcomes, report, opts, "Dummies", "Dummy variables")


def run_presence(doc, opts):
    s = opts.settings
    threshold = float(s.get("threshold", 0.0))
    empty_absent = _flag(s, "empty_absent")
    inp = read_input(doc, opts)
    outcomes = _run_each(inp, lambda name, v: _coded(
        name, coding.presence_absence(v, threshold, empty_absent), "{name} present"),
        numeric=False)
    report = _report("Presence/absence coding", inp, [], outcomes)
    return _write(doc, inp, outcomes, report, opts, "Presence-absence", "Presence/absence")


def run_bins(doc, opts):
    s = opts.settings
    method, k = s.get("method", "width"), int(s.get("k", 4))
    funcs = {"width": coding.bins_equal_width, "quantile": coding.bins_quantile,
             "kmeans": coding.bins_kmeans}
    inp = read_input(doc, opts)
    outcomes = _run_each(inp, lambda name, v: _coded(name, funcs[method](v, k), "{name} bin"))
    report = _report("Discretisation into bins", inp,
                     [("Bins are", "right-closed: (from, to]; the first bin includes its lower "
                                   "limit")] if method != "kmeans" else [], outcomes)
    return _write(doc, inp, outcomes, report, opts, "Bins", "Bins")


def run_recode(doc, opts):
    rules = coding.parse_rules(str(opts.settings.get("rules", "")))
    inp = read_input(doc, opts)
    outcomes = _run_each(inp, lambda name, v: _coded(name, coding.recode(v, rules),
                                                      "{name} recoded"), numeric=False)
    report = _report("Recoded categories", inp, [], outcomes)
    return _write(doc, inp, outcomes, report, opts, "Recoded", "Recode categories")


def run_ranks(doc, opts):
    s = opts.settings
    want_ranks, want_scores = _flag(s, "ranks", True), _flag(s, "scores")
    ties, method = s.get("ties", "average"), s.get("score_method", "blom")
    if not (want_ranks or want_scores):
        raise calc.RangeError("Choose ranks, normal scores, or both.")

    def func(name, values):
        cols, params = [], []
        if want_ranks:
            r = transform.ranks(values, ties)
            cols.append((f"rank({name})", r.values))
            params += r.params
        if want_scores:
            ns = transform.normal_scores(values, method)
            cols.append((f"nscore({name})", ns.values))
            params += ns.params
        return cols, params, [], []

    inp = read_input(doc, opts)
    outcomes = _run_each(inp, func)
    report = _report("Ranks and normal scores", inp,
                     [("Ranks", "1 = smallest value")] if want_ranks else [], outcomes)
    return _write(doc, inp, outcomes, report, opts, "Ranks", "Ranks and normal scores")


RUNNERS = {"transform": run_transform, "scale": run_scale, "logratio": run_logratio,
           "dummies": run_dummies, "presence": run_presence, "bins": run_bins,
           "recode": run_recode, "ranks": run_ranks}
