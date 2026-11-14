"""Analyse a screening design sheet made by LibreStats and write the report.

The design table is recognised by its headers: "Point type", "Block", the
"<factor> (coded)" columns, the actual-level columns that follow them, and the
response columns after those. The report holds the effects table, Lenth's
method, the curvature check, a Pareto chart of the standardized effects and a
half-normal plot of the effects (both native Calc charts).
"""

import math
from dataclasses import dataclass

import uno
from com.sun.star.awt import Point, Size

from . import calc
from . import influence
from .distributions import t_isf
from . import screening as sc
from .viz_calc import LINE_SOLID, _add_chart, _axis_titles, _series, \
    _write_block, insert_scatter_chart

CODED = " (coded)"
P_FORMAT = '[<0.0001]"<0.0001";0.0000'
G = "0.####"
BAR_COLOR = 0x004586
LINE_COLORS = (0xC9211E, 0xFF8000)       # ME / critical t, SME
PARETO_SIZE = (14000, 9000)
HALF_NORMAL_SIZE = (14000, 10000)
REFERENCES = ["Lenth, R. V. (1989). Technometrics, 31(4), 469-473.",
              "Daniel, C. (1959). Technometrics, 1(4), 311-341.",
              "Montgomery, D. C. (2017). Design and Analysis of Experiments (9th ed.), "
              "chapters 6-8. Wiley."]
COOK_REFERENCES = ["Cook, R. D. (1977). Technometrics, 19(1), 15-18.",
                   "Bollen, K. A., & Jackman, R. W. (1990). In J. Fox & J. S. Long (Eds.), "
                   "Modern Methods of Data Analysis. Sage.",
                   "Hoaglin, D. C., & Welsch, R. E. (1978). The American Statistician, 32(1), "
                   "17-22.",
                   "Weisberg, S. (2014). Applied Linear Regression (4th ed.). Wiley."]
LEVEL_FORMAT = "0.0%"
INFLUENCE_SIZE = (14000, 10000)     # 1/100 mm
INFLUENCE_MIN, INFLUENCE_RANGE = 120, 650   # marker diameter, 1/100 mm
INFLUENCE_LABELS = 3                         # label the runs with the largest D and |t*|
WORTH_A_LOOK = 2.0                           # influence plot: |t*| reference line
REF_GREY = 0x808080
OUTLIER_ALPHA = 0.05                         # Bonferroni line: two-sided, over all runs


def add_diagnostics(add, table, observations, diags, p, df, fmt_y):
    """Write the per-run diagnostics (influence.py) through the report's add(): twelve
    columns (A-L). Returns the row index of the table header and the number of runs,
    for influence_plot()."""
    n = len(diags)
    cutoff = influence.leverage_cutoff(n, p)
    add("Diagnostics per run: outliers, leverage and influence", bold=True, overflow=True)
    add(f"r = residual / (S sqrt(1 - h)) (internally studentized); t* = r sqrt((n - p - 1) / "
        f"(n - p - r^2)) (externally studentized, t({df - 1}) for a run that is not an "
        f"outlier); p Bonferroni = min(1, n p) for testing all {n} runs; Cook's D = "
        f"r^2 h / (p (1 - h)), p = {p} model terms; F level = F({p}, {df}) probability at D.",
        overflow=True)
    lev = (f"high leverage: h > 2p/n = {cutoff:.3g}" if cutoff is not None else
           "high leverage: not assessed (p/n > 0.5)")
    add("Flags: outlier = Bonferroni p < 0.05; check = p < 0.05 only; " + lev +
        "; influential = F level 50% or more (D at or above the median of F).",
        overflow=True)
    header = add("Std order" if table.std else "Row", "Observed", "Fitted", "Residual",
                 "Stud. resid. r", "Leverage h", "Ext. stud. t*", "p (outlier)",
                 "p Bonferroni", "Cook's D", "F level", "Flags", bold=True)
    blank = lambda v: "" if v is None else v  # noqa: E731
    lines, no_cooks, no_external = [], 0, 0
    for i, d in zip(observations, diags):
        label = table.std[i] if table.std else i + 2
        no_cooks += d.cooks is None
        no_external += d.external is None and d.cooks is not None
        add(label, d.observed, d.fitted, d.residual, blank(d.studentized), d.leverage,
            blank(d.external), blank(d.p_outlier), blank(d.p_bonferroni), blank(d.cooks),
            blank(d.level), ", ".join(d.flags()) if d.cooks is not None else "n/a",
            fmt=(None, fmt_y, fmt_y, fmt_y, "0.000", "0.000", "0.000", P_FORMAT, P_FORMAT,
                 "0.0000", LEVEL_FORMAT))
        line = influence.verdict(f"{label:g}" if isinstance(label, float) else label, d, n, p)
        if line:
            lines.append(line)
    add("Flagged runs (by " + ("standard order" if table.std else "row") + ")" +
        ("" if lines else ": none"), bold=True, overflow=True)
    for line in lines:
        add(line, overflow=True)
    if cutoff is None:
        add(f"Warning: p/n = {p}/{n} > 0.5: the model has nearly as many terms as runs, "
            "so leverage is high everywhere and the 2p/n rule cannot pick out unusual runs.",
            overflow=True)
    if df - 1 <= 0:
        add("Note: the outlier test needs n - p - 1 > 0 (at least 2 residual degrees of "
            "freedom); t* and its p values are n/a.", overflow=True)
    elif no_external:
        add(f"Note: t* is n/a for {no_external} run(s): without that run the other runs "
            "fit exactly (r^2 = n - p), so its t* is infinite.", overflow=True)
    if no_cooks:
        add(f"Note: the diagnostics are n/a for {no_cooks} run(s) with leverage 1 (the model "
            "cannot be fitted without them) or when the fit is perfect (S = 0).",
            overflow=True)
    shown = {i for i, d in enumerate(diags) if d.flags()}
    for key in (lambda i: diags[i].cooks or 0.0, lambda i: abs(diags[i].external or 0.0)):
        shown |= set(sorted(range(n), key=key, reverse=True)[:INFLUENCE_LABELS])
    return header, n, cutoff, sorted(shown), p, df


def influence_plot(doc, sheet, header, n, cutoff, labelled, p, df, model_name="",
                   size=INFLUENCE_SIZE, greyscale=False):
    """Native XY chart of the diagnostics table: leverage (x) against t* (y), marker
    area proportional to Cook's D. The flagged runs and those with the largest D and
    |t*| are labelled with their run (design points often share a leverage, so labels
    on every point would pile up). Dashed reference lines: |t*| = 2 ("worth a look"),
    +- the Bonferroni critical t (two-sided 0.05 over all n runs, n - p - 1 df) and
    h = 2p/n; their end points are written in a small table (column N) and the chart
    is placed right of it, below any chart or picture already in that area."""
    data = sheet.getCellRangeByPosition(5, header, 6, header + n)       # h, t*
    labels = [sheet.getCellByPosition(0, header + 1 + i).getString() for i in range(n)]
    cooks = [sheet.getCellByPosition(9, header + 1 + i) for i in range(n)]
    cooks = [c.getValue() if c.getType().value == "VALUE" else None for c in cooks]
    if not any(c is not None for c in cooks):
        return None
    ts = [sheet.getCellByPosition(6, header + 1 + i) for i in range(n)]
    ts = [abs(c.getValue()) for c in ts if c.getType().value == "VALUE"]
    crit = t_isf(OUTLIER_ALPHA / (2 * n), df - 1) if df - 1 > 0 else None
    top = max(ts + [WORTH_A_LOOK] + ([crit] if crit else [])) * 1.1
    top = math.ceil(top * 2) / 2                        # a round axis limit

    # reference lines: (name, colour, (x0, y0), (x1, y1), in the legend)
    w = WORTH_A_LOOK
    grey_c, crit_c, lev_c = (REF_GREY, 0x000000, 0x555555) if greyscale else \
        (REF_GREY, LINE_COLORS[0], LINE_COLORS[1])
    lines = [(f"|t*| = {w:g} (worth a look)", grey_c, (0, w), (1, w), True),
             (f"t* = -{w:g}", grey_c, (0, -w), (1, -w), False)]
    if crit:
        lines += [(f"Bonferroni critical t = \u00b1{crit:.2f}", crit_c, (0, crit),
                   (1, crit), True),
                  (f"t* = -{crit:.2f}", crit_c, (0, -crit), (1, -crit), False)]
    if cutoff is not None:
        lines.append((f"2p/n = {cutoff:.3g} (high leverage)", lev_c, (cutoff, -top),
                      (cutoff, top), True))
    col = 13
    table = [["Influence plot reference lines", "x", "y"]]
    spans = []
    for name, _, a, b, _ in lines:
        spans.append((header + len(table), header + len(table) + 1))
        table += [[name, *a], [name, *b]]
    _write_block(sheet, col, header, table)
    sheet.getCellRangeByPosition(col, header, col + 2, header).CharWeight = calc.BOLD

    pos = sheet.getCellByPosition(col + 4, header).Position
    left = pos.X
    page = sheet.getDrawPage()
    for j in range(page.getCount()):
        shape = page.getByIndex(j)
        sp, sz = shape.Position, shape.Size
        if sp.X + sz.Width > left and sp.X < left + size[0]:
            pos = Point(pos.X, max(pos.Y, sp.Y + sz.Height + 500))
    title = "Influence plot" + (f": {model_name}" if model_name else "")
    chart = insert_scatter_chart(doc, sheet, "Influence plot", data, pos, size,
                                 title + " (marker area ~ Cook's D)", "Leverage h",
                                 "Externally studentized residual t*",
                                 0x333333 if greyscale else BAR_COLOR)
    diagram = chart.getDiagram()
    x_axis, y_axis = diagram.getXAxis(), diagram.getYAxis()
    x_axis.CrossoverPosition = uno.Enum("com.sun.star.chart.ChartAxisPosition", "START")
    x_axis.AutoMin = x_axis.AutoMax = False
    x_axis.Min, x_axis.Max = 0.0, 1.0
    y_axis.AutoMin = y_axis.AutoMax = False
    y_axis.Min, y_axis.Max = -top, top
    diagram.HasYAxisGrid = False           # the reference lines take over that role
    chart_type = chart.getFirstDiagram().getCoordinateSystems()[0].getChartTypes()[0]
    series = chart_type.getDataSeries()[0]
    series.ShowLegendEntry = False
    biggest = max(c for c in cooks if c is not None) or 1.0
    for i, c in enumerate(cooks):
        point = series.getDataPointByIndex(i)
        symbol = point.Symbol
        symbol.Style = uno.Enum("com.sun.star.chart2.SymbolStyle", "STANDARD")
        symbol.StandardSymbol = 8          # circle (LibreOffice: 0 square, 1 diamond, 8 circle)
        d = INFLUENCE_MIN + INFLUENCE_RANGE * math.sqrt(max(c or 0.0, 0.0) / biggest)
        symbol.Size = Size(int(d), int(d))
        point.Symbol = symbol
    _label_points(chart, [(i, labels[i]) for i in labelled])
    for i in labelled:                                  # labels beside the markers
        series.getDataPointByIndex(i).LabelPlacement = uno.getConstantByName(
            "com.sun.star.chart.DataLabelPlacement.RIGHT")
    for (_, color, _, _, legend), (first, last) in zip(lines, spans):
        _add_line_series(chart, chart_type, sheet, col, first, last, color, legend)
    chart.HasLegend = True
    chart.getLegend().Alignment = uno.Enum("com.sun.star.chart.ChartLegendPosition", "BOTTOM")
    _plain_axes(chart)
    return chart


def _add_line_series(chart, chart_type, sheet, col, first, last, color, legend=True):
    """A dashed line through the (x, y) rows first..last of the reference table (name in
    column col, x and y in the next two)."""
    ctx = uno.getComponentContext()
    smgr = ctx.ServiceManager
    provider = chart.getDataProvider()

    def labelled(c, role):
        rng = sheet.getCellRangeByPosition(c, first, c, last).AbsoluteName
        values = provider.createDataSequenceByRangeRepresentation(rng)
        values.Role = role
        seq = smgr.createInstanceWithContext("com.sun.star.chart2.data.LabeledDataSequence",
                                             ctx)
        seq.setValues(values)
        seq.setLabel(provider.createDataSequenceByRangeRepresentation(
            sheet.getCellByPosition(col, first).AbsoluteName))
        return seq

    s = smgr.createInstanceWithContext("com.sun.star.chart2.DataSeries", ctx)
    s.setData((labelled(col + 1, "values-x"), labelled(col + 2, "values-y")))
    chart_type.addDataSeries(s)
    s.Color = color
    s.LineWidth = 35
    s.LineStyle = uno.Enum("com.sun.star.drawing.LineStyle", "DASH")
    dash = uno.createUnoStruct("com.sun.star.drawing.LineDash")
    dash.Style = uno.Enum("com.sun.star.drawing.DashStyle", "RECT")
    dash.Dots, dash.DotLen, dash.Dashes, dash.DashLen, dash.Distance = 0, 0, 1, 200, 150
    s.LineDash = dash
    symbol = s.Symbol
    symbol.Style = uno.Enum("com.sun.star.chart2.SymbolStyle", "NONE")
    s.Symbol = symbol
    s.ShowLegendEntry = legend


@dataclass
class ScreeningOptions:
    range_text: str = ""          # the design table (header row included)
    has_labels: bool = True       # always true for a design sheet
    response: str = ""            # empty: the first response column
    model: str = "2fi"            # 2fi | main
    alpha: float = 0.05
    charts: bool = True


@dataclass
class DesignTable:
    factors: list                 # factor names
    coded: list                   # rows of coded values
    ptype: list
    blocks: list
    responses: dict               # name -> list of values (None = empty)
    actual: list = None           # rows of actual levels (None if any is text)
    std: list = None              # standard order of each row (None: not numeric)
    labels: list = None           # headers of the actual-level columns, e.g. "Time (min)"


def read_design_table(data, screening=True):
    """Parse a LibreStats design table (data array with the header row first).
    screening=True: two-level designs only (coded -1, 0, +1); False: any design, e.g.
    a response surface design with axial points."""
    if not data:
        raise sc.ScreeningError("The design table is empty.")
    header = [str(h).strip() for h in data[0]]
    end = header.index("") if "" in header else len(header)
    header = header[:end]
    coded_cols = [j for j, h in enumerate(header) if h.endswith(CODED)]
    if "Point type" not in header or "Block" not in header or len(coded_cols) < 2:
        raise sc.ScreeningError("Select a design sheet made with DOE / RSM > New Design: it "
                                "needs the columns Point type, Block and the '(coded)' factor "
                                "columns.")
    k = len(coded_cols)
    first_response = coded_cols[-1] + 1 + k
    names = header[first_response:]
    if not names:
        raise sc.ScreeningError("The design table has no response columns.")
    pt, bl = header.index("Point type"), header.index("Block")
    so = header.index("Std order") if "Std order" in header else None
    rows = [r for r in data[1:] if any(str(v).strip() for v in r[:end])]

    def number(v, what, row):
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return v
        raise sc.ScreeningError(f"Row {row + 2}: {what} must be a number.")

    # point types first: a response surface design fails the coded-value check below
    # (axial points at +-alpha), and that message would not tell what is wrong
    ptype = [int(number(r[pt], "Point type", i)) for i, r in enumerate(rows)]
    if screening and any(p in sc.RSM_POINT_TYPES for p in ptype):
        raise sc.ScreeningError(sc.RSM_DESIGN)
    coded, blocks, actual = [], [], []
    for i, r in enumerate(rows):
        values = [number(r[j], header[j], i) for j in coded_cols]
        if screening:
            bad = [header[j] for j, v in zip(coded_cols, values) if v not in (-1, 0, 1)]
            if bad:
                raise sc.ScreeningError(f"Row {i + 2}: {bad[0]} must be -1, 0 or +1.")
            values = [int(v) for v in values]
        coded.append(values)
        blocks.append(number(r[bl], "Block", i))
        levels = [r[coded_cols[-1] + 1 + j] for j in range(k)]
        if actual is not None and all(isinstance(v, (int, float)) and not isinstance(v, bool)
                                      for v in levels):
            actual.append([float(v) for v in levels])
        else:
            actual = None
    responses = {}
    for offset, name in enumerate(names):
        col = first_response + offset
        values = []
        for i, r in enumerate(rows):
            v = r[col]
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                values.append(float(v))
            elif str(v).strip():
                raise sc.ScreeningError(f"{name}, row {i + 2}: '{v}' is not a number.")
            else:
                values.append(None)
        responses[name] = values
    factors = [header[j][:-len(CODED)] for j in coded_cols]
    std = [r[so] for r in rows] if so is not None else None
    if std is not None and not all(isinstance(v, (int, float)) and not isinstance(v, bool)
                                   for v in std):
        std = None
    labels = [header[coded_cols[-1] + 1 + j] for j in range(k)]
    return DesignTable(factors, coded, ptype, blocks, responses, actual, std, labels)


def default_table_range(doc):
    """The design table on the active sheet (region around A1), or ''."""
    sheet = doc.getCurrentController().getActiveSheet()
    if sheet.getCellByPosition(0, 0).getString() != "Std order":
        return ""
    cursor = sheet.createCursorByRange(sheet.getCellByPosition(0, 0))
    cursor.collapseToCurrentRegion()
    return cursor.AbsoluteName


def response_names(doc, range_text, screening=True):
    """(response names, factor names, problem): the names, or [], [] and why the table
    cannot be read."""
    try:
        table = read_design_table(calc.used_data(calc.resolve_range(doc, range_text)),
                                  screening)
        return list(table.responses), table.factors, None
    except (calc.RangeError, sc.ScreeningError) as exc:
        return [], [], str(exc)


def _fmt(x, digits=4):
    return "" if x is None else f"{x:.{digits}g}"


def _term_label(term, factors):
    letters = sc.effect_name(term.mask)
    if len(letters) == 1:
        return f"{letters} ({factors[sc.LETTERS.index(letters)]})"
    return letters


def _factors_text(term, factors):
    return " x ".join(factors[i] for i in range(len(factors)) if term.mask >> i & 1)


def run_screening(doc, opts):
    cell_range = calc.resolve_range(doc, opts.range_text)
    table = read_design_table(calc.used_data(cell_range))
    response = opts.response.strip() or next(iter(table.responses))
    if response not in table.responses:
        raise sc.ScreeningError(f"There is no response column '{response}'. Choose one of: "
                                + ", ".join(table.responses) + ".")
    a = sc.analyse(table.coded, table.ptype, table.blocks, table.responses[response],
                   opts.model, opts.alpha)
    factors = table.factors

    fitted, long = [], []          # (row index, values, bold, formats); long rows overflow
    scale = max([abs(t.effect) for t in a.terms] + [abs(e) for _, e in a.contrasts] + [1e-12])
    E = "0." + "0" * min(8, max(2, 4 - math.floor(math.log10(scale))))   # effects, in units
    T, P, R = "0.000", P_FORMAT, "0.0000"

    def add(*values, bold=False, overflow=False, fmt=()):
        index = len(fitted) + len(long)
        (long if overflow else fitted).append((index, list(values), bold, fmt))
        return index

    add(f"Screening analysis of {response}", bold=True, overflow=True)
    add("Design table", cell_range.AbsoluteName, overflow=True)
    add("Runs analysed", a.n)
    if a.missing:
        add("Runs without a response", a.missing)
    add("Model", "main effects and two-factor interactions" if opts.model == "2fi"
        else "main effects", overflow=True)
    if len(set(table.blocks)) > 1:
        add("Blocks", len(set(table.blocks)))
    error = {"residual": f"residual mean square ({a.df} df)",
             "lenth": f"Lenth's pseudo standard error (d = m/3 = {a.lenth.df:.3g} df)"
             if a.lenth else "", "none": "none (saturated model)"}[a.error]
    add("Error estimate", error, overflow=True)
    if a.mse is not None:
        add("S (residual standard deviation)", a.mse ** 0.5, fmt=(None, E))
    if a.r2 is not None:
        add("R-squared", a.r2, fmt=(None, R))
    if a.r2_adj is not None:
        add("Adjusted R-squared", a.r2_adj, fmt=(None, R))
    add()
    add("Effects", bold=True)
    add("Term", "Factor(s)", "Effect", "Coefficient", "SE of effect",
        "t (Lenth)" if a.error == "lenth" else "t", "p", "Significant", bold=True)
    add("Constant", "", "", a.intercept, fmt=(None, None, None, E))
    crit = a.critical_t(opts.alpha)
    for t in a.terms:
        sig = "" if t.p is None else "yes" if t.p < opts.alpha else "no"
        add(t.name, _factors_text(t, factors), t.effect, t.coef, t.se, t.t, t.p, sig,
            fmt=(None, None, E, E, E, T, P))
    aliased = [t for t in a.terms if t.aliases]
    if aliased:
        add()
        add("Aliases (each effect above estimates the sum of its alias chain)", bold=True,
            overflow=True)
        for t in aliased:
            add(f"{t.name} {('+ ' if not t.alias_text().startswith('-') else '')}"
                f"{t.alias_text()}", overflow=True)
    if a.block_effects:
        add()
        add("Confounded with blocks: " + ", ".join(t.name + (" " + t.alias_text()
                                                             if t.aliases else "")
                                                    for t in a.block_effects), overflow=True)
    if a.lenth:
        L = a.lenth
        add()
        add("Lenth's method", bold=True, overflow=True)
        add("Contrasts (m)", L.m)
        add("Pseudo standard error (PSE)", L.pse, fmt=(None, E))
        add(f"Margin of error ME (alpha = {opts.alpha:g})", L.me, fmt=(None, E))
        add("Simultaneous margin of error SME", L.sme, fmt=(None, E))
        beyond_me = [t.name for t, e in a.contrasts if abs(e) > L.me]
        beyond_sme = [t.name for t, e in a.contrasts if abs(e) > L.sme]
        add("Effects beyond ME", ", ".join(beyond_me) or "none", overflow=True)
        add("Effects beyond SME", ", ".join(beyond_sme) or "none", overflow=True)
    if a.pure_df and a.df > a.pure_df and a.mse is not None:
        lof_df = a.df - a.pure_df
        lof_ss = a.mse * a.df - a.pure_ss
        f = (lof_ss / lof_df) / (a.pure_ss / a.pure_df) if a.pure_ss > 0 else None
        add()
        add("Lack of fit", bold=True, overflow=True)
        add("Pure error (repeated runs): SS, df", a.pure_ss, a.pure_df, fmt=(None, G))
        add("Lack of fit: SS, df", lof_ss, lof_df, fmt=(None, G))
        if f is not None:
            add("F, p", f, sc.f_sf(f, lof_df, a.pure_df), fmt=(None, T, P))
    c = a.curvature
    if c:
        add()
        add("Curvature check (centre points)", bold=True, overflow=True)
        add("Mean of factorial runs", c.mean_factorial, fmt=(None, E))
        add("Mean of centre points", c.mean_centre, fmt=(None, E))
        add("Centre minus factorial (block-adjusted)", c.estimate, fmt=(None, E))
        add("SS curvature (1 df)", c.ss, fmt=(None, G))
        if c.f is not None:
            add("F, p", c.f, c.p, fmt=(None, T, P))
            verdict = ("Significant curvature: the response is not linear in the factors; a "
                       "response surface design is needed to model it."
                       if c.p < opts.alpha else
                       "No significant curvature: a first-order model (with interactions) "
                       "describes the region.")
            add(verdict, overflow=True)
    add()
    diag = None
    if a.diagnostics is not None:
        diag = add_diagnostics(add, table, a.observations, a.diagnostics, a.n - a.df, a.df,
                               E)
    else:
        add("Diagnostics per run (Cook's distance): not available, because the model has "
            "no residual degrees of freedom.", overflow=True)
    if a.notes:
        add()
        for n in a.notes:
            add("Note: " + n, overflow=True)
    add()
    add("References", bold=True, overflow=True)
    for r in REFERENCES + COOK_REFERENCES:
        add(r, overflow=True)

    with calc.undo_context(doc, "LibreStats: Screening analysis"):
        sheet = calc.new_sheet(doc, "Screening analysis")
        formats = _Formats(doc)
        for index, values, bold, fmt in fitted:
            if values:
                _write_block(sheet, 0, index, [values])
            for j, code in enumerate(fmt):
                if code:
                    sheet.getCellByPosition(j, index).NumberFormat = formats.key(code)
            if bold:
                sheet.getCellRangeByPosition(0, index, 11, index).CharWeight = calc.BOLD
        sheet.getCellRangeByPosition(0, 0, 11, 0).getColumns().OptimalWidth = True
        for index, values, bold, _ in long:
            if values:
                _write_block(sheet, 0, index, [values])
            if bold:
                sheet.getCellByPosition(0, index).CharWeight = calc.BOLD
        if opts.charts:
            _charts(doc, sheet, a, factors, opts.alpha, crit, formats)
            if diag:
                influence_plot(doc, sheet, *diag,
                               model_name="main effects and two-factor interactions"
                               if opts.model == "2fi" else "main effects")
    doc.getCurrentController().setActiveSheet(sheet)
    return sheet.getName()


class _Formats:
    def __init__(self, doc):
        self.formats = doc.getNumberFormats()
        self.locale = uno.createUnoStruct("com.sun.star.lang.Locale")

    def key(self, code):
        key = self.formats.queryKey(code, self.locale, False)
        return key if key != -1 else self.formats.addNew(code, self.locale)


def _charts(doc, sheet, a, factors, alpha, crit, formats):
    col = 13                                   # chart tables start in column N
                                               # (the run diagnostics use A-L)
    row = 0
    terms = [t for t in a.terms if t.t is not None]
    origin = sheet.getCellByPosition(col + 5, 0).Position
    if terms:
        terms.sort(key=lambda t: abs(t.t))      # smallest first: drawn at the bottom
        lines = []
        if a.error == "lenth":
            lines = [("ME", a.lenth.me / a.lenth.pse), ("SME", a.lenth.sme / a.lenth.pse)]
        elif crit is not None:
            lines = [(f"t crit ({alpha:g})", crit)]
        header = ["Term", "|t|"] + [name for name, _ in lines]
        rows = [header] + [[_term_label(t, factors), abs(t.t)] + [v for _, v in lines]
                           for t in terms]
        table = _write_block(sheet, col, row, rows)
        sheet.getCellRangeByPosition(col, row, col + len(header) - 1, row).CharWeight = \
            calc.BOLD
        sheet.getCellRangeByPosition(col + 1, row + 1, col + len(header) - 1,
                                     row + len(terms)).NumberFormat = formats.key("0.000")
        what = "Lenth's pseudo t" if a.error == "lenth" else "t"
        refs = ", ".join(f"{n} = {v:.3g}" for n, v in lines)
        _pareto(doc, sheet, table, origin, f"Pareto chart of standardized effects",
                f"|{what}|" + (f"   (lines: {refs})" if refs else ""), len(lines))
        row += len(rows) + 1

    if a.lenth and a.contrasts:
        L = a.lenth
        ordered = sorted(a.contrasts, key=lambda c: abs(c[1]))
        q = sc.half_normal_quantiles(len(ordered))
        rows = [["Half-normal quantile", "|Effect|", "Lenth line (PSE x q)", "Effect"]]
        rows += [[qi, abs(e), L.pse * qi, t.name] for qi, (t, e) in zip(q, ordered)]
        table = _write_block(sheet, col, row, rows)
        sheet.getCellRangeByPosition(col, row, col + 3, row).CharWeight = calc.BOLD
        sheet.getCellRangeByPosition(col, row + 1, col + 2,
                                     row + len(ordered)).NumberFormat = formats.key("0.000")
        data = sheet.getCellRangeByPosition(col, row, col + 2, row + len(ordered))
        pos = Point(origin.X, origin.Y + PARETO_SIZE[1] + 500)
        chart = insert_scatter_chart(doc, sheet, "Half-normal plot", data, pos,
                                     HALF_NORMAL_SIZE, "Half-normal plot of the effects",
                                     "Half-normal quantile", "|Effect|", BAR_COLOR,
                                     LINE_COLORS[0])
        labels = [(i, t.name) for i, (t, e) in enumerate(ordered) if abs(e) > L.me]
        _label_points(chart, labels)
        _plain_axes(chart)
    sheet.getCellRangeByPosition(col, 0, col + 3, 0).getColumns().OptimalWidth = True


def _pareto(doc, sheet, table, pos, title, value_title, n_lines):
    """Horizontal bars of |t|, largest on top, with reference lines."""
    chart = _add_chart(doc, sheet, "Pareto chart", table, pos, PARETO_SIZE, True, title,
                       "", value_title)
    _axis_titles(chart, "Term", value_title)
    cs = chart.getFirstDiagram().getCoordinateSystems()[0]
    bars, series = _series(chart)
    series = list(series)
    bars.setDataSeries(tuple(series[:1]))
    bars.GapwidthSequence = (60,)
    series[0].Color = BAR_COLOR
    if n_lines:
        ctx = uno.getComponentContext()
        line = ctx.ServiceManager.createInstanceWithContext(
            "com.sun.star.chart2.LineChartType", ctx)
        line.setDataSeries(tuple(series[1:]))
        for s, color in zip(series[1:], LINE_COLORS):
            s.Color = color
            s.LineWidth = 50
            s.LineStyle = LINE_SOLID
            symbol = s.Symbol
            symbol.Style = uno.Enum("com.sun.star.chart2.SymbolStyle", "NONE")
            s.Symbol = symbol
        cs.addChartType(line)
    cs.SwapXAndYAxis = True
    _plain_axes(chart)
    chart.getDiagram().getYAxisTitle().TextRotation = 0     # the value axis is now horizontal
    return chart


def _plain_axes(chart):
    """Axis labels in the chart's own number format, not that of the source cells."""
    diagram = chart.getDiagram()
    for axis in (diagram.getXAxis(), diagram.getYAxis()):
        axis.LinkNumberFormatToSource = False


def _label_points(chart, labels):
    """Text labels on points of the first series."""
    if not labels:
        return
    ctx = uno.getComponentContext()
    ct = chart.getFirstDiagram().getCoordinateSystems()[0].getChartTypes()[0]
    s = ct.getDataSeries()[0]
    for i, text in labels:
        point = s.getDataPointByIndex(i)
        f = ctx.ServiceManager.createInstanceWithContext(
            "com.sun.star.chart2.DataPointCustomLabelField", ctx)
        f.FieldType = uno.Enum("com.sun.star.chart2.DataPointCustomLabelFieldType", "TEXT")
        f.setString(text)
        point.CustomLabelFields = (f,)
        label = point.Label
        label.ShowCustomLabel = True
        label.ShowNumber = False
        point.Label = label


RUNNERS = {"doeAnalyseScreening": run_screening}
