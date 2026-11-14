"""Analyse a response surface design sheet and write the report.

The design table is read as for the screening analysis (section 9.3 of the
documentation) but may hold any coded values, e.g. axial points at +-alpha.
The report holds the coefficients, the ANOVA, the equation in actual units,
the stationary point with its canonical analysis, the best and worst predicted
response within the experimental region, a response map (a grid of predicted
values coloured by height) and the fitted surface as a 3D picture.
"""

import math
import re
from dataclasses import dataclass

from . import calc
from . import reduction
from . import rsm
from . import rsm_plots
from . import surface as sf
from .doe import LETTERS
from .screening_calc import (COOK_REFERENCES, G, P_FORMAT, _Formats, add_diagnostics,
                             influence_plot, read_design_table)
from com.sun.star.awt import Point

from .viz_calc import _write_block, insert_svg

MAP_SIZE = 21                    # grid points per axis of the response map
MAP_FORMAT_DIGITS = 4
SURFACE_WIDTH = 15000            # 1/100 mm
PX = 22                          # 1/100 mm per SVG pixel of the new plots
MAP_COLUMN, MAP_STEP = 10, 23    # first response map in column K; maps side by side
VIF_HIGH = 10                    # VIF above this: flagged (Montgomery, 2017)
MODEL_NAMES = {"quadratic": "full quadratic", "interaction": "linear + two-factor interactions",
               "linear": "linear", "cubic": "reduced cubic (full cubic without aliased terms)"}
REFERENCES = ["Box, G. E. P., & Wilson, K. B. (1951). Journal of the Royal Statistical Society "
              "B, 13(1), 1-45.",
              "Myers, R. H., Montgomery, D. C., & Anderson-Cook, C. M. (2016). Response "
              "Surface Methodology (4th ed.). Wiley.",
              "Montgomery, D. C. (2017). Design and Analysis of Experiments (9th ed.), "
              "chapter 11. Wiley.",
              "Hurvich, C. M., & Tsai, C.-L. (1989). Biometrika, 76(2), 297-307.",
              "Draper, N. R., & Smith, H. (1998). Applied Regression Analysis (3rd ed.). "
              "Wiley."]


@dataclass
class RSMOptions:
    range_text: str = ""          # the design table (header row included)
    has_labels: bool = True
    response: str = ""            # empty: the first response column
    model: str = "quadratic"      # quadratic | interaction | linear
    alpha: float = 0.05
    x_factor: str = ""            # one response map of these two factors (automation); empty:
    y_factor: str = ""            # every factor pair, at most rsm.MAX_MAPS (largest range)
    charts: bool = True
    reduction: str = "none"       # none | manual | backward | subsets (reduction.METHODS)
    reduction_alpha: float = 0.05     # backward elimination: alpha-out
    hierarchy: bool = True
    locked: tuple = ()            # term names never removed
    remove: tuple = ()            # manual: term names to remove
    top: int = 5                  # all subsets: models shown
    rank: int = 1                 # all subsets: the model used
    exclude_aliased: bool = False     # reduced cubic: leave aliased stand-ins out
    goal: str = "max"             # max | min: the optimum used by the plots
    fixed: str = "optimum"        # maps and surfaces: factors not shown at optimum |
                                  # centre | custom
    fixed_values: str = ""        # custom: "Time=13; EhOH=50" (names or letters, real units)
    residuals: str = "raw"        # residuals-by-factor plot: raw | studentized
    anchor: str = "optimum"       # optimum trace: optimum | centre (perturbation plot)
    greyscale: bool = False       # journal style: grey shades, black lines


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
    return _monomial_label(rsm.monomial(term), factors)


def _monomial_label(m, factors):
    """(0, 0, 1) -> 'Temperature^2 x Time'."""
    return " x ".join(factors[i] + (f"^{m.count(i)}" if m.count(i) > 1 else "")
                      for i in sorted(set(m)))


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
    try:
        red = reduction.reduce(table.coded, table.blocks, table.responses[response],
                               opts.model, opts.reduction, opts.reduction_alpha,
                               opts.hierarchy, opts.locked, opts.remove, opts.top, opts.rank,
                               opts.exclude_aliased)
    except reduction.ReductionError as exc:
        raise rsm.RSMError(str(exc)) from None
    fit = red.final
    marks = rsm.alias_marks(fit.aliases)
    term_label = lambda name: f"{name} ({marks[name]})" if name in marks else name  # noqa
    blocked = any(t.kind == "block" for t in fit.terms)
    if opts.x_factor.strip() or opts.y_factor.strip():
        ix = _factor_index(opts.x_factor, factors, 0)
        iy = _factor_index(opts.y_factor, factors, 1 if ix != 1 else 0)
        pairs, every = [(ix, iy)], None
    else:
        pairs, every = rsm.map_pairs(fit)
        ix, iy = pairs[0]
    if ix == iy:
        raise rsm.RSMError("Choose two different factors for the response map.")

    def actual(j, coded):
        return centres[j] + halves[j] * coded

    labels = table.labels or factors
    model_text = MODEL_NAMES[opts.model].split(" (")[0] + \
        (", reduced" if opts.reduction != "none" else "")
    terms_plain = ", ".join(t.name for t in fit.terms if t.kind not in ("constant", "block"))
    subtitle = f"{model_text}: {terms_plain or 'constant only'}"
    goal_best = fit.maximum if opts.goal == "max" else fit.minimum
    goal_name = "maximum" if opts.goal == "max" else "minimum"
    fixed, fixed_name, plot_notes = _fixed_point(opts, goal_best, goal_name, factors,
                                                 centres, halves, k)
    ys = [v for v in table.responses[response] if v is not None]
    scale = max(abs(v) for v in ys) or 1.0
    Y = "0." + "0" * min(8, max(2, MAP_FORMAT_DIGITS - math.floor(math.log10(scale))))
    C, T, P, R = "0.0000", "0.000", P_FORMAT, "0.0000"
    EQ = "0.0#######"          # actual-unit coefficients span many magnitudes
    fitted, long = [], []

    def add(*values, bold=False, overflow=False, fmt=()):
        index = len(fitted) + len(long)
        (long if overflow else fitted).append((index, list(values), bold, fmt))
        return index

    add(f"Response surface analysis of {response}", bold=True, overflow=True)
    add("Design table", cell_range.AbsoluteName, overflow=True)
    add("Runs analysed", fit.n)
    if fit.missing:
        add("Runs without a response", fit.missing)
    add("Model", MODEL_NAMES[opts.model] + (" (coded units)" if opts.reduction == "none" else
                                            ", reduced (coded units)"), overflow=True)
    terms_text = ", ".join(term_label(t.name) for t in fit.terms
                           if t.kind not in ("constant", "block"))
    add("Terms in the model", (terms_text or "constant only") + (
        " (+ blocks)" if any(t.kind == "block" for t in fit.terms) else ""), overflow=True)
    add("S (residual standard deviation)", fit.mse ** 0.5, fmt=(None, Y))
    if fit.r2 is not None:
        add("R-squared", fit.r2, fmt=(None, R))
        add("Adjusted R-squared", fit.r2_adj, fmt=(None, R))
    if fit.r2_pred is not None:
        add("Predicted R-squared", fit.r2_pred, fmt=(None, R))
    if fit.r2_max is not None:
        add("R-squared ceiling (pure error)", fit.r2_max, fmt=(None, R))
    if opts.reduction != "none":
        _reduction_section(add, red, MODEL_NAMES[opts.model], Y, term_label)
    add()
    add("Coefficients (coded units)", bold=True, overflow=True)
    add("Term", "Factor(s)", "Coefficient", "SE", "t", "p", "Significant", "VIF", bold=True)
    for t in fit.terms:
        add(term_label(t.name), _term_factors(t, factors), t.coef, t.se, t.t, t.p,
            "" if t.kind == "constant" else "yes" if t.p < opts.alpha else "no",
            "" if t.vif is None else t.vif, fmt=(None, None, C, C, T, P, None, "0.00"))
    high = [t.name for t in fit.terms if t.vif is not None and t.vif > VIF_HIGH]
    if high:
        add(f"VIF above {VIF_HIGH}: {', '.join(high)}. These terms are strongly correlated "
            "with other terms in this design, so their coefficients are imprecise and can "
            "change a lot when terms are added or removed.", overflow=True)
    if opts.model == "cubic":
        add()
        add("Reduced cubic: aliased terms removed", bold=True, overflow=True)
        n_blocks = sum(t.kind == "block" for t in fit.terms)
        incl = f" (including {n_blocks} block term{'s' if n_blocks > 1 else ''})" \
            if n_blocks else ""
        supported = fit.full_terms - len(fit.aliases)
        add(f"The design supports {supported} of the {fit.full_terms} full-cubic terms{incl}"
            f" ({len(fit.aliases)} aliased term(s) removed); the final model has "
            f"{len(fit.terms)}. A removed term cannot be told apart from the terms it is "
            "aliased with: its effect is carried by their coefficients, so a term marked "
            "'aliased with' estimates the sum of both. (Coded units; terms are kept in "
            "hierarchical order: linear, square, interaction, then ABC, A^2B and A^3 types.)",
            overflow=True)
        if red.excluded:
            add("Excluded from the model as aliased (option): " + ", ".join(red.excluded)
                + ".", overflow=True)
        for a in fit.aliases:
            add(a.text(), overflow=True)
        if not fit.aliases:
            add("No term is aliased: the full cubic can be estimated.", overflow=True)
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
    add("Equation in actual units" + (" (average block)" if blocked else ""), bold=True,
        overflow=True)
    equation = rsm.actual_equation(fit, centres, halves)
    add("Term", "Factor(s)", "Coefficient", bold=True)
    for name, coef in equation:
        add(name, _monomial_label(rsm.name_monomial(name), factors), coef, fmt=(None, None, EQ))
    add()
    add("Factor levels: actual = centre + half-range x coded", bold=True, overflow=True)
    add("Factor", "Name", "Centre (0)", "Half-range", "Coded range in the design", bold=True)
    for j in range(k):
        lo, hi = fit.bounds[j]
        add(LETTERS[j], factors[j], centres[j], halves[j], f"{lo:.4g} to {hi:.4g}",
            fmt=(None, None, G, G))

    observed = [d.observed for d in fit.diagnostics]
    unrun = rsm.corners_not_run(fit.bounds, [table.coded[i] for i in fit.observations])
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
            beyond = _beyond_observed(st.response, observed, "predicted response at the "
                                      "stationary point")
            if beyond:
                add(beyond, overflow=True)
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
        add(f"{label} predicted response within the experimental region"
            + (" (found numerically)" if opts.model == "cubic" else ""), bold=True,
            overflow=True)
        add("Factor", "Name", "Coded", "Actual", bold=True)
        for j in range(k):
            add(LETTERS[j], factors[j], best[0][j], actual(j, best[0][j]),
                fmt=(None, None, C, G))
        add("Predicted response", "", best[1], fmt=(None, None, Y))
        beyond = _beyond_observed(best[1], observed, f"predicted {label.lower()}")
        if beyond:
            add(beyond, overflow=True)
        if any(all(abs(a - b) <= 1e-6 for a, b in zip(best[0], c)) for c in unrun):
            add(f"Warning: the {label.lower()} lies at a corner of the region that was not "
                "run: an extrapolation. Confirm it with a run there.", overflow=True)
    add()
    diag = add_diagnostics(add, table, fit.observations, fit.diagnostics, fit.n - fit.df,
                           fit.df, Y)
    add()
    if opts.charts:
        names = lambda prs: ", ".join(f"{factors[a]} x {factors[b]}" for a, b in prs)  # noqa
        text = f"Response maps and 3D surfaces (to the right): {names(pairs)}"
        if every and len(every) > len(pairs):
            text += (f"; the {len(pairs)} of the {len(every)} factor pairs with the largest "
                     f"range of predicted values. Not drawn: "
                     f"{names([p for p in every if p not in pairs])}")
        add(text + ".", overflow=True)
        add("Plots (to the right, below the maps): observed vs predicted, residuals by "
            "factor (" + ("externally studentized" if opts.residuals == "studentized"
                          else "raw") + "), optimum trace (anchor: " +
            (f"the predicted {goal_name}" if opts.anchor == "optimum" and goal_best else
             "the centre") + "), influence plot" +
            (" (greyscale)" if opts.greyscale else "") + ".", overflow=True)
        for line in plot_notes:
            add(line, overflow=True)
    corners = 2 ** k
    if not unrun:
        add("The experimental region is the box of the coded settings that were run. All "
            f"{corners} corners of the box were run, so predictions inside it interpolate "
            "between measured settings; outside it they are extrapolations.", overflow=True)
    else:
        add("The experimental region is the box of the coded settings that were run. "
            f"{len(unrun)} of its {corners} corners were not run (e.g. axial points reach "
            "beyond the factorial cube, or a Box-Behnken design has no corners), so "
            "predictions near those corners are extrapolations.", overflow=True)
    for n in fit.notes:
        add("Note: " + n, overflow=True)
    add()
    add("References", bold=True, overflow=True)
    for r in REFERENCES + COOK_REFERENCES:
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
                sheet.getCellRangeByPosition(0, index, 11, index).CharWeight = calc.BOLD
        sheet.getCellRangeByPosition(0, 0, 11, 0).getColumns().OptimalWidth = True
        for index, values, bold, _ in long:
            if values:
                _write_block(sheet, 0, index, [values])
            if bold:
                sheet.getCellByPosition(0, index).CharWeight = calc.BOLD
        if opts.charts:
            for j, (a, b) in enumerate(pairs):
                _response_map(doc, sheet, fit, factors, a, b, actual, response, formats, Y,
                              col=MAP_COLUMN + j * MAP_STEP, fixed=fixed,
                              fixed_name=fixed_name, labels=labels, subtitle=subtitle,
                              greyscale=opts.greyscale)
            _model_plots(doc, sheet, fit, table, opts, response, labels, centres, halves,
                         actual, subtitle, goal_best, goal_name)
            influence_plot(doc, sheet, *diag, model_name=model_text, greyscale=opts.greyscale)
    doc.getCurrentController().setActiveSheet(sheet)
    return sheet.getName()


def _beyond_observed(value, observed, what):
    """A warning when a predicted value lies outside the measured responses, else None."""
    lo, hi = min(observed), max(observed)
    if value < lo:
        return (f"Warning: the {what} ({value:.4g}) is below the lowest measured response "
                f"({lo:.4g}): the model extrapolates beyond the data; confirm with runs.")
    if value > hi:
        return (f"Warning: the {what} ({value:.4g}) is above the highest measured response "
                f"({hi:.4g}): the model extrapolates beyond the data; confirm with runs.")
    return None


def _na(v):
    return "n/a" if v is None else v


def _reduction_section(add, red, model_name, fmt_y, term_label=lambda n: n):
    """Model reduction: settings, steps or the top-N models, start vs final model."""
    R, P = "0.0000", P_FORMAT
    add()
    head = f"Model reduction: {reduction.METHOD_NAMES[red.method]}"
    if red.method == "backward":
        head += f", alpha = {red.alpha:g}"
    head += ", hierarchy " + ("kept" if red.hierarchy else "off")
    add(head, bold=True, overflow=True)
    if red.locked:
        add("Locked terms (never removed)", ", ".join(red.locked), overflow=True)
    if red.method == "backward":
        if red.steps:
            add("Step", "Term removed", "p at removal", bold=True)
            for i, (t, p) in enumerate(red.steps, 1):
                add(i, term_label(t), p, fmt=(None, None, P))
        else:
            add(f"No removable term had p > {red.alpha:g}: nothing was removed.",
                overflow=True)
    elif red.method == "manual":
        add("Removed (chosen)", ", ".join(term_label(t) for t in red.removed) or "none",
            overflow=True)
    elif red.method == "subsets":
        add(f"{red.evaluated} models fitted"
            + (" (hierarchical)" if red.hierarchy else "")
            + "; ranked by AICc = n ln(SSE/n) + 2p + 2p(p+1)/(n-p-1), lowest first (ties: "
            "fewer terms, then higher predicted R-squared).", overflow=True)
        add("Rank", "AICc", "Delta AICc", "R-squared", "Adj. R-sq.", "Pred. R-sq.",
            "Lack-of-fit p", "Used", "Terms", bold=True)
        best = red.ranked[0].aicc
        shown = list(red.candidates)
        if red.chosen not in shown:
            shown.append(red.chosen)
        for c in shown:
            rank = next(i for i, x in enumerate(_all_ranked(red), 1) if x is c)
            add(rank, c.aicc, c.aicc - best, c.r2, c.r2_adj, _na(c.r2_pred), _na(c.lof_p),
                "<-" if c is red.chosen else "",
                ", ".join(term_label(t) for t in c.terms) or "constant only",
                fmt=(None, "0.000", "0.000", R, R, R, P))
        add("Delta AICc = AICc - lowest AICc; models with Delta AICc < 2 are about equally "
            "supported by the data (Burnham & Anderson, 2002).", overflow=True)
        if red.rank != 1:
            add(f"The model of rank {red.rank} is used (chosen in the dialog).",
                overflow=True)
    add("", "Starting model", "Final model", bold=True)
    s, f = red.start, red.final
    for label, a, b, code in (("R-squared", s.r2, f.r2, R),
                              ("Adjusted R-squared", s.r2_adj, f.r2_adj, R),
                              ("Predicted R-squared", s.r2_pred, f.r2_pred, R),
                              ("Lack-of-fit p", s.lof_p, f.lof_p, P),
                              ("AICc", s.aicc, f.aicc, "0.000"),
                              ("S", s.mse ** 0.5, f.mse ** 0.5, fmt_y),
                              ("Terms (with the constant)", len(s.terms), len(f.terms), None)):
        add(label, _na(a), _na(b), fmt=(None, code, code))
    add(red.summary(model_name), overflow=True)
    add("Note: p-values of a selected model are optimistic; judge it by predicted R-squared "
        "and confirmation runs.", overflow=True)
    for w in red.warnings:
        add("Warning: " + w, overflow=True)
    removed = ", ".join(red.removed)
    add("To re-run identically: model reduction 'Manual', removing " +
        (removed if removed else "nothing") + ".", overflow=True)


def _all_ranked(red):
    """The candidates in rank order as far as they were ranked (top N plus the chosen)."""
    return red.ranked


def _response_map(doc, sheet, fit, factors, ix, iy, actual, response, formats, fmt,
                  col=10, row=0, fixed=None, fixed_name="the centre", labels=None,
                  subtitle="", greyscale=False):
    """Grid of predicted values over factors ix (columns) and iy (rows), the others at
    their centre, coloured by height; and the fitted surface as a 3D picture below it."""
    k = len(factors)
    (x0, x1), (y0, y1) = fit.bounds[ix], fit.bounds[iy]
    xs = [x0 + (x1 - x0) * i / (MAP_SIZE - 1) for i in range(MAP_SIZE)]
    ys = [y1 - (y1 - y0) * i / (MAP_SIZE - 1) for i in range(MAP_SIZE)]   # high at top

    def point(cx, cy):
        x = list(fixed) if fixed is not None else [0.0] * k
        x[ix], x[iy] = cx, cy
        return x

    grid = [[fit.predict(point(cx, cy)) for cx in xs] for cy in ys]
    lo = min(min(r) for r in grid)
    hi = max(max(r) for r in grid)
    labels = labels or factors
    fixed = fixed if fixed is not None else [0.0] * k
    others = [f"{labels[j]} fixed at {fixed_name} = {actual(j, fixed[j]):.4g}"
              for j in range(k) if j not in (ix, iy)]
    footnote = "; ".join(others)
    title = [[f"Response map: predicted {response}" + (f" ({footnote})" if footnote else "")]]
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
            cell.CellBackColor = int((sf.grey(t) if greyscale else sf.color(t))[1:], 16)
            cell.CharColor = 0xFFFFFF if t < 0.6 else 0x000000
    block.getColumns().OptimalWidth = True
    # the long title (and the note below) are written after the widths are set, so they
    # overflow to the right instead of widening the first column of the map
    _write_block(sheet, col, row, title)
    sheet.getCellByPosition(col, row).CharWeight = calc.BOLD
    note_row = row + 3 + MAP_SIZE
    _write_block(sheet, col, note_row,
                 [["Colours: viridis scale from the lowest (dark blue) to the highest "
                   "(yellow) predicted value on the map."]])

    rows = [(actual(ix, cx), actual(iy, cy), v) for cy, r in zip(ys, grid)
            for cx, v in zip(xs, r)]
    surface = sf.build_surface(rows)
    svg = sf.render_svg(surface, title=f"Fitted {response}", x_title=labels[ix],
                        y_title=labels[iy], z_title=response, greyscale=greyscale,
                        subtitle=_short(subtitle, 100), footnote=footnote)
    pos = sheet.getCellByPosition(col, note_row + 2).Position
    insert_svg(doc, sheet, f"Fitted surface: {factors[ix]} x {factors[iy]}", svg, pos,
               SURFACE_WIDTH)


def _short(text, n):
    text = rsm_plots.pretty(text)
    return text if len(text) <= n else text[:n - 3] + "..."


def _fixed_point(opts, best, goal_name, factors, centres, halves, k):
    """Coded values for the factors not shown on a map or surface, a description, and
    report notes."""
    notes = []
    if opts.fixed == "optimum":
        if best is not None:
            return list(best[0]), f"the predicted {goal_name}", notes
        notes.append(f"Note: no predicted {goal_name} was found; the factors not shown are "
                     "held at their centre.")
        return [0.0] * k, "the centre", notes
    if opts.fixed == "centre":
        return [0.0] * k, "the centre", notes
    point = [0.0] * k
    given = set()
    for part in opts.fixed_values.replace(",", ";").split(";"):
        if not part.strip():
            continue
        if "=" not in part:
            raise rsm.RSMError(f"Fixed values: write name=value, e.g. 'Time=13' (got "
                               f"'{part.strip()}').")
        name, value = (v.strip() for v in part.split("=", 1))
        j = _factor_index(name, factors, None)
        try:
            real = float(value.replace(",", "."))
        except ValueError:
            raise rsm.RSMError(f"Fixed values: '{value}' is not a number.") from None
        point[j] = (real - centres[j]) / halves[j]
        given.add(j)
    missing = [factors[j] for j in range(k) if j not in given]
    if missing:
        notes.append("Note: no fixed value given for " + ", ".join(missing) +
                     "; held at the centre on the maps and surfaces.")
    return point, "chosen values", notes


def _model_plots(doc, sheet, fit, table, opts, response, labels, centres, halves, actual,
                 subtitle, best, goal_name):
    """Observed vs predicted, residuals by factor and the optimum trace, as pictures in a
    row below the response maps and surfaces."""
    page = sheet.getDrawPage()
    bottom = max((page.getByIndex(i).Position.Y + page.getByIndex(i).Size.Height
                  for i in range(page.getCount())), default=0)
    x = sheet.getCellByPosition(MAP_COLUMN, 0).Position.X
    y = bottom + 800
    diags = fit.diagnostics
    std = [table.std[i] if table.std else i + 2 for i in fit.observations]
    marks = [f"{s:g}" if isinstance(s, float) else str(s) for s in std]

    def place(name, svg):
        nonlocal x
        width = int(PX * float(re.search(r'width="([\d.]+)"', svg).group(1)))
        insert_svg(doc, sheet, name, svg, Point(x, y), width)
        x += width + 600

    flagged = [m if d.outlier or d.influential else None for m, d in zip(marks, diags)]
    place("Observed vs predicted", rsm_plots.observed_vs_predicted(
        [d.observed for d in diags], [d.fitted for d in diags], flagged, fit.r2,
        fit.r2_pred, f"Observed vs predicted {response}", subtitle, response, opts.greyscale))

    segments = []
    for j in range(fit.k):
        levels = sorted({round(table.coded[i][j], 9) for i in fit.observations})
        index = [levels.index(round(table.coded[i][j], 9)) for i in fit.observations]
        # tick text: the real setting as written on the Design sheet (no rounding noise)
        real = {}
        for i in fit.observations:
            real.setdefault(round(table.coded[i][j], 9), table.actual[i][j])
        segments.append((labels[j], [(lv, f"{real[lv]:.4g}") for lv in levels], index))
    studentized = opts.residuals == "studentized"
    values = [d.external if studentized else d.residual for d in diags]
    place("Residuals by factor", rsm_plots.residuals_by_factor(
        segments, values, f"Residuals of {response} by factor", subtitle,
        "Externally studentized residual t*" if studentized else f"Residual ({response})",
        opts.greyscale, guides=studentized))

    anchor = list(best[0]) if opts.anchor == "optimum" and best is not None else [0.0] * fit.k
    anchor_name = f"predicted {goal_name}" if opts.anchor == "optimum" and best else "centre"
    svg, _, _ = rsm_plots.optimum_trace(
        fit, anchor, labels, centres, halves,
        ("Optimum trace" if anchor_name != "centre" else "Perturbation plot") + f": {response}",
        subtitle, response, anchor_name, opts.greyscale)
    place("Optimum trace", svg)


RUNNERS = {"doeAnalyseRSM": run_rsm}
