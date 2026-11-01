"""Analyse a response surface design sheet and write the report.

The design table is read as for the screening analysis (section 9.3 of the
documentation) but may hold any coded values, e.g. axial points at +-alpha.
The report holds the coefficients, the ANOVA, the equation in actual units,
the stationary point with its canonical analysis, the best and worst predicted
response within the experimental region, a response map (a grid of predicted
values coloured by height) and the fitted surface as a 3D picture.
"""

import math
from dataclasses import dataclass

from . import calc
from . import rsm
from . import surface as sf
from .doe import LETTERS
from .screening_calc import G, P_FORMAT, _Formats, read_design_table
from .viz_calc import _write_block, insert_svg

MAP_SIZE = 21                    # grid points per axis of the response map
MAP_FORMAT_DIGITS = 4
SURFACE_WIDTH = 15000            # 1/100 mm
MODEL_NAMES = {"quadratic": "full quadratic", "interaction": "linear + two-factor interactions",
               "linear": "linear"}
REFERENCES = ["Box, G. E. P., & Wilson, K. B. (1951). Journal of the Royal Statistical Society "
              "B, 13(1), 1-45.",
              "Myers, R. H., Montgomery, D. C., & Anderson-Cook, C. M. (2016). Response "
              "Surface Methodology (4th ed.). Wiley.",
              "Montgomery, D. C. (2017). Design and Analysis of Experiments (9th ed.), "
              "chapter 11. Wiley."]


@dataclass
class RSMOptions:
    range_text: str = ""          # the design table (header row included)
    has_labels: bool = True
    response: str = ""            # empty: the first response column
    model: str = "quadratic"      # quadratic | interaction | linear
    alpha: float = 0.05
    x_factor: str = ""            # response map axes: factor names; empty = first two
    y_factor: str = ""
    charts: bool = True


def _factor_index(name, factors, default):
    if not name.strip():
        return default
    if name in factors:
        return factors.index(name)
    if name.upper() in LETTERS[:len(factors)]:
        return LETTERS.index(name.upper())
    raise rsm.RSMError(f"There is no factor '{name}'. Choose one of: " + ", ".join(factors)
                       + ".")


def _term_factors(term, factors):
    if term.kind in ("constant", "block"):
        return ""
    if term.kind == "square":
        return f"{factors[term.factors[0]]}^2"
    return " x ".join(factors[i] for i in term.factors)


def run_rsm(doc, opts):
    cell_range = calc.resolve_range(doc, opts.range_text)
    table = read_design_table(calc.used_data(cell_range), screening=False)
    response = opts.response.strip() or next(iter(table.responses))
    if response not in table.responses:
        raise rsm.RSMError(f"There is no response column '{response}'. Choose one of: "
                           + ", ".join(table.responses) + ".")
    if table.actual is None:
        raise rsm.RSMError("Response surface analysis needs numeric actual levels for every "
                           "factor (the columns after the coded ones).")
    factors, k = table.factors, len(table.factors)
    centres, halves, misfit = rsm.coding(table.coded, table.actual)
    if misfit > 0.01:
        raise rsm.RSMError("The actual levels do not match the coded levels (actual = centre "
                           "+ half-range x coded). Check the factor columns.")
    fit = rsm.fit_model(table.coded, table.blocks, table.responses[response], opts.model)
    ix = _factor_index(opts.x_factor, factors, 0)
    iy = _factor_index(opts.y_factor, factors, 1 if ix != 1 else 0)
    if ix == iy:
        raise rsm.RSMError("Choose two different factors for the response map.")

    def actual(j, coded):
        return centres[j] + halves[j] * coded

    ys = [v for v in table.responses[response] if v is not None]
    scale = max(abs(v) for v in ys) or 1.0
    Y = "0." + "0" * min(8, max(2, MAP_FORMAT_DIGITS - math.floor(math.log10(scale))))
    C, T, P, R = "0.0000", "0.000", P_FORMAT, "0.0000"
    EQ = "0.0#######"          # actual-unit coefficients span many magnitudes
    fitted, long = [], []

    def add(*values, bold=False, overflow=False, fmt=()):
        (long if overflow else fitted).append((len(fitted) + len(long), list(values), bold,
                                               fmt))

    add(f"Response surface analysis of {response}", bold=True, overflow=True)
    add("Design table", cell_range.AbsoluteName, overflow=True)
    add("Runs analysed", fit.n)
    if fit.missing:
        add("Runs without a response", fit.missing)
    add("Model", MODEL_NAMES[opts.model] + " (coded units)", overflow=True)
    add("S (residual standard deviation)", fit.mse ** 0.5, fmt=(None, Y))
    if fit.r2 is not None:
        add("R-squared", fit.r2, fmt=(None, R))
        add("Adjusted R-squared", fit.r2_adj, fmt=(None, R))
    if fit.r2_pred is not None:
        add("Predicted R-squared", fit.r2_pred, fmt=(None, R))
    add()
    add("Coefficients (coded units)", bold=True, overflow=True)
    add("Term", "Factor(s)", "Coefficient", "SE", "t", "p", "Significant", bold=True)
    for t in fit.terms:
        add(t.name, _term_factors(t, factors), t.coef, t.se, t.t, t.p,
            "" if t.kind == "constant" else "yes" if t.p < opts.alpha else "no",
            fmt=(None, None, C, C, T, P))
    add()
    add("Analysis of variance", bold=True, overflow=True)
    add("Source", "DF", "Adj SS", "Adj MS", "F", "p", bold=True)
    for r in fit.anova:
        ms = r.ss / r.df if r.df else None
        add(r.source, r.df, r.ss, ms if r.source != "Total" else "",
            "" if r.f is None else r.f, "" if r.p is None else r.p,
            fmt=(None, None, G, G, T, P))
    lof = next((r for r in fit.anova if r.source.strip() == "Lack of fit"), None)
    if lof and lof.p is not None:
        add("Lack of fit is significant: the model does not describe the response well "
            "enough; consider a transformation or more terms." if lof.p < opts.alpha else
            "Lack of fit is not significant: the model fits within the pure error.",
            overflow=True)
    add()
    add(f"Equation in actual units (average block)", bold=True, overflow=True)
    equation = rsm.actual_equation(fit, centres, halves)
    add("Term", "Factor(s)", "Coefficient", bold=True)
    for name, coef in equation:
        idx = [LETTERS.index(c) for c in name.replace("^2", "") if c in LETTERS] \
            if name != "Constant" else []
        label = "" if not idx else factors[idx[0]] + "^2" if name.endswith("^2") else \
            " x ".join(factors[i] for i in idx)
        add(name, label, coef, fmt=(None, None, EQ))
    add()
    add("Factor levels: actual = centre + half-range x coded", bold=True, overflow=True)
    add("Factor", "Name", "Centre (0)", "Half-range", "Coded range in the design", bold=True)
    for j in range(k):
        lo, hi = fit.bounds[j]
        add(LETTERS[j], factors[j], centres[j], halves[j], f"{lo:.4g} to {hi:.4g}",
            fmt=(None, None, G, G))

    st = fit.stationary
    if st is not None:
        add()
        add("Stationary point and canonical analysis", bold=True, overflow=True)
        if st.coded is not None:
            add("Factor", "Name", "Coded", "Actual", bold=True)
            for j in range(k):
                add(LETTERS[j], factors[j], st.coded[j], actual(j, st.coded[j]),
                    fmt=(None, None, C, G))
            add("Predicted response", "", st.response, fmt=(None, None, Y))
            add("Nature", st.nature, overflow=True)
            add("Inside the experimental region", "yes" if st.inside else "no")
        add("Eigenvalues", *st.eigenvalues, fmt=(None,) + (C,) * k)
        verdict = {
            "maximum": "All eigenvalues are negative: the stationary point is a maximum.",
            "minimum": "All eigenvalues are positive: the stationary point is a minimum.",
            "saddle point": "The eigenvalues have mixed signs: the stationary point is a "
                            "saddle point (neither maximum nor minimum); the optimum lies on "
                            "the boundary of the region (see below).",
        }.get(st.nature, "")
        if st.ridge:
            verdict += (" An eigenvalue is close to 0: the surface is a ridge, and settings "
                        "along its direction give nearly the same response.")
        if st.coded is not None and not st.inside:
            verdict += (" The stationary point lies outside the experimental region; the "
                        "model is an extrapolation there.")
        add(verdict.strip(), overflow=True)
    for label, best in (("Maximum", fit.maximum), ("Minimum", fit.minimum)):
        if best is None:
            continue
        add()
        add(f"{label} predicted response within the experimental region", bold=True,
            overflow=True)
        add("Factor", "Name", "Coded", "Actual", bold=True)
        for j in range(k):
            add(LETTERS[j], factors[j], best[0][j], actual(j, best[0][j]),
                fmt=(None, None, C, G))
        add("Predicted response", "", best[1], fmt=(None, None, Y))
    add()
    add("The experimental region is the box of the coded settings that were run; for a "
        "rotatable design its corners lie beyond the axial points, so predictions there are "
        "extrapolations.", overflow=True)
    for n in fit.notes:
        add("Note: " + n, overflow=True)
    add()
    add("References", bold=True, overflow=True)
    for r in REFERENCES:
        add(r, overflow=True)

    with calc.undo_context(doc, "LibreStats: Response surface analysis"):
        sheet = calc.new_sheet(doc, "RSM analysis")
        formats = _Formats(doc)
        for index, values, bold, fmt in fitted:
            if values:
                _write_block(sheet, 0, index, [values])
            for j, code in enumerate(fmt):
                if code and j < len(values) and isinstance(values[j], (int, float)):
                    sheet.getCellByPosition(j, index).NumberFormat = formats.key(code)
            if bold:
                sheet.getCellRangeByPosition(0, index, 7, index).CharWeight = calc.BOLD
        sheet.getCellRangeByPosition(0, 0, 7, 0).getColumns().OptimalWidth = True
        for index, values, bold, _ in long:
            if values:
                _write_block(sheet, 0, index, [values])
            if bold:
                sheet.getCellByPosition(0, index).CharWeight = calc.BOLD
        if opts.charts:
            _response_map(doc, sheet, fit, factors, ix, iy, actual, response, formats, Y)
    doc.getCurrentController().setActiveSheet(sheet)
    return sheet.getName()


def _response_map(doc, sheet, fit, factors, ix, iy, actual, response, formats, fmt):
    """Grid of predicted values over factors ix (columns) and iy (rows), the others at
    their centre, coloured by height; and the fitted surface as a 3D picture."""
    col, row = 10, 0
    k = len(factors)
    (x0, x1), (y0, y1) = fit.bounds[ix], fit.bounds[iy]
    xs = [x0 + (x1 - x0) * i / (MAP_SIZE - 1) for i in range(MAP_SIZE)]
    ys = [y1 - (y1 - y0) * i / (MAP_SIZE - 1) for i in range(MAP_SIZE)]   # high at top

    def point(cx, cy):
        x = [0.0] * k
        x[ix], x[iy] = cx, cy
        return x

    grid = [[fit.predict(point(cx, cy)) for cx in xs] for cy in ys]
    lo = min(min(r) for r in grid)
    hi = max(max(r) for r in grid)
    others = [f"{factors[j]} = {actual(j, 0):.4g}" for j in range(k) if j not in (ix, iy)]
    title = [[f"Response map: predicted {response}" +
              (f" (other factors at their centre: {', '.join(others)})" if others else "")]]
    _write_block(sheet, col, row, title)
    sheet.getCellByPosition(col, row).CharWeight = calc.BOLD
    header = [f"{factors[iy]} \\ {factors[ix]}"] + [actual(ix, cx) for cx in xs]
    body = [[actual(iy, cy)] + r for cy, r in zip(ys, grid)]
    block = _write_block(sheet, col, row + 1, [header] + body)
    sheet.getCellRangeByPosition(col, row + 1, col + MAP_SIZE, row + 1).CharWeight = calc.BOLD
    sheet.getCellRangeByPosition(col, row + 2, col, row + 1 + MAP_SIZE).CharWeight = calc.BOLD
    sheet.getCellRangeByPosition(col + 1, row + 1, col + MAP_SIZE,
                                 row + 1).NumberFormat = formats.key(G)
    sheet.getCellRangeByPosition(col, row + 2, col, row + 1 + MAP_SIZE).NumberFormat = \
        formats.key(G)
    cells = sheet.getCellRangeByPosition(col + 1, row + 2, col + MAP_SIZE, row + 1 + MAP_SIZE)
    cells.NumberFormat = formats.key(fmt)
    for r, values in enumerate(grid):
        for c, v in enumerate(values):
            t = 0.5 if hi == lo else (v - lo) / (hi - lo)
            cell = cells.getCellByPosition(c, r)
            cell.CellBackColor = int(sf.color(t)[1:], 16)
            cell.CharColor = 0xFFFFFF if t < 0.6 else 0x000000
    block.getColumns().OptimalWidth = True
    note_row = row + 3 + MAP_SIZE
    _write_block(sheet, col, note_row,
                 [["Colours: viridis scale from the lowest (dark blue) to the highest "
                   "(yellow) predicted value on the map."]])

    rows = [(actual(ix, cx), actual(iy, cy), v) for cy, r in zip(ys, grid)
            for cx, v in zip(xs, r)]
    surface = sf.build_surface(rows)
    svg = sf.render_svg(surface, title=f"Fitted {response}", x_title=factors[ix],
                        y_title=factors[iy], z_title=response)
    pos = sheet.getCellByPosition(col, note_row + 2).Position
    insert_svg(doc, sheet, "Fitted surface", svg, pos, SURFACE_WIDTH)


RUNNERS = {"doeAnalyseRSM": run_rsm}
