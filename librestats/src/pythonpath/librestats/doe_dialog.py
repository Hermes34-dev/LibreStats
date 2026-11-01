"""Dialogs for new screening and response surface designs (DOE / RSM menu)."""

from . import calc
from . import dialogs as d
from .doe_calc import DesignOptions
from .rsm_calc import RSMOptions
from .screening_calc import ScreeningOptions, default_table_range, response_names

KIND_RADIOS = (("full", "kindFull"), ("fractional", "kindFractional"), ("pb", "kindPB"))


def _whole(text, what, low):
    value = d.number(text, what)
    if value is None or value != int(value) or value < low:
        raise ValueError(f"{what} must be a whole number of {low} or more.")
    return int(value)


def show_design(ctx, doc, parent, kind="full"):
    """Show the dialog; returns DesignOptions, or None if cancelled."""
    title = "LibreStats - New Screening Design"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 290)
    b.label("lblRange", 10, 8, 55, "Factor table:")
    b.edit("range", 70, 8, 180, range_text)
    b.check("labels", 10, 24, 240, "First row contains column headings", labels)
    b.label("lblHint", 10, 38, 240, "Columns: Name, Unit, Low, High (Unit may be left out). "
            "Levels may be text.")

    b.line("lineDesign", 54, "Design")
    b.radio("kindFull", 10, 66, 240, "Full factorial (2 levels): all low/high combinations",
            kind == "full")
    b.radio("kindFractional", 10, 80, 240, "Fractional factorial 2^(k-p), minimum aberration",
            kind == "fractional")
    b.radio("kindPB", 10, 94, 240, "Plackett-Burman (8, 12, 16, 20 or 24 runs)", kind == "pb")
    b.label("lblRuns", 10, 110, 85, "Runs:")
    b.edit("runs", 100, 110, 40, "")
    b.label("lblRunsHint", 145, 110, 105, "empty = automatic")

    b.line("lineOptions", 128, "Options")
    b.label("lblCentre", 10, 140, 85, "Centre points per block:")
    b.edit("centre", 100, 140, 40, "0")
    b.label("lblReplicates", 10, 156, 85, "Replicates:")
    b.edit("replicates", 100, 156, 40, "1")
    b.label("lblBlocks", 10, 172, 85, "Blocks:")
    b.edit("blocks", 100, 172, 40, "1")
    b.label("lblBlocksHint", 145, 172, 105, "replicates x 1, 2, 4 or 8")
    b.check("randomize", 10, 188, 120, "Randomise the run order", True)
    b.label("lblSeed", 130, 188, 40, "Seed:")
    b.edit("seed", 170, 188, 50, "")
    b.label("lblSeedHint", 130, 202, 120, "(empty = new random seed)")

    b.line("lineResponses", 218, "Responses")
    b.label("lblResponses", 10, 230, 85, "Response columns:")
    b.edit("responses", 100, 230, 150, "Response")
    b.label("lblResponsesHint", 10, 244, 240, "Several names separated by commas, e.g. "
            "Yield, Purity.")
    b.ok_cancel(266)

    dialog = b.create(parent)
    problems = {}

    def read(dlg):
        text = lambda name: dlg.getControl(name).getText()  # noqa: E731
        opts = DesignOptions(
            range_text=text("range"), has_labels=d.checked(dlg, "labels"),
            kind=next(k for k, name in KIND_RADIOS if d.checked(dlg, name)),
            randomize=d.checked(dlg, "randomize"), responses=text("responses"))
        problems.pop("msg", None)
        try:
            if text("runs").strip():
                opts.runs = _whole(text("runs"), "The number of runs", 2)
            opts.centre_points = _whole(text("centre") or "0", "Centre points", 0)
            opts.replicates = _whole(text("replicates") or "1", "Replicates", 1)
            opts.blocks = _whole(text("blocks") or "1", "Blocks", 1)
            if text("seed").strip():
                opts.seed = _whole(text("seed"), "The random seed", 0)
        except ValueError as exc:
            problems["msg"] = str(exc)
        return opts

    return d.run(ctx, doc, parent, dialog, title, read, lambda _: problems.get("msg"))


ALPHA_RADIOS = (("rotatable", "alphaRotatable"), ("orthogonal", "alphaOrthogonal"),
                ("face", "alphaFace"), ("custom", "alphaCustom"))


def show_rsm(ctx, doc, parent, kind="ccd"):
    """Show the dialog; returns DesignOptions, or None if cancelled."""
    title = "LibreStats - New Response Surface Design"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 352)
    b.label("lblRange", 10, 8, 55, "Factor table:")
    b.edit("range", 70, 8, 180, range_text)
    b.check("labels", 10, 24, 240, "First row contains column headings", labels)
    b.label("lblHint", 10, 38, 240, "Columns: Name, Unit, Low, High (Unit may be left out). "
            "Numeric levels.")

    b.line("lineDesign", 54, "Design")
    b.radio("kindCCD", 10, 66, 240, "Central composite (CCD): factorial + axial + centre "
            "points", kind == "ccd")
    b.radio("kindBB", 10, 80, 240, "Box-Behnken (3 to 7 factors): 3 levels, no extreme "
            "corners", kind == "bbd")

    b.line("lineCCD", 96, "Central composite only")
    b.label("lblAlpha", 10, 108, 85, "Axial distance alpha:")
    b.radio("alphaRotatable", 100, 108, 150, "rotatable", True)
    b.radio("alphaOrthogonal", 100, 122, 150, "orthogonal blocking (needs blocks)", False)
    b.radio("alphaFace", 100, 136, 150, "face-centred (alpha = 1)", False)
    b.radio("alphaCustom", 100, 150, 50, "custom:", False)
    b.edit("customAlpha", 150, 150, 40, "")
    b.check("fractionalCube", 10, 166, 240, "Resolution V fraction as factorial part "
            "(5 or more factors)", True)
    b.check("inscribed", 10, 180, 240, "Inscribed: axial points at Low and High (limits "
            "that cannot be exceeded)", False)

    b.line("lineOptions", 196, "Options (leave centre points empty for the usual numbers)")
    b.label("lblCentre", 10, 208, 120, "Centre points (per cube block):")
    b.edit("centre", 135, 208, 40, "")
    b.label("lblAxialCentre", 10, 224, 120, "Centre points in the axial block:")
    b.edit("axialCentre", 135, 224, 40, "")
    b.label("lblReplicates", 10, 240, 120, "Replicates:")
    b.edit("replicates", 135, 240, 40, "1")
    b.label("lblBlocks", 10, 256, 120, "Blocks:")
    b.edit("blocks", 135, 256, 40, "1")
    b.label("lblBlocksHint", 10, 270, 240, "Replicates x 1, 2, 3 or 5 (CCD); x 1, or 3 for "
            "a 4-factor Box-Behnken.")
    b.check("randomize", 10, 286, 120, "Randomise the run order", True)
    b.label("lblSeed", 130, 286, 40, "Seed:")
    b.edit("seed", 170, 286, 50, "")

    b.line("lineResponses", 302, "Responses")
    b.label("lblResponses", 10, 314, 85, "Response columns:")
    b.edit("responses", 100, 314, 150, "Response")
    b.ok_cancel(334)

    dialog = b.create(parent)
    problems = {}

    def read(dlg):
        text = lambda name: dlg.getControl(name).getText()  # noqa: E731
        opts = DesignOptions(
            range_text=text("range"), has_labels=d.checked(dlg, "labels"),
            kind="ccd" if d.checked(dlg, "kindCCD") else "bbd",
            alpha=next(a for a, name in ALPHA_RADIOS if d.checked(dlg, name)),
            fractional_cube=d.checked(dlg, "fractionalCube"),
            inscribed=d.checked(dlg, "inscribed"), centre_points=None,
            randomize=d.checked(dlg, "randomize"), responses=text("responses"))
        problems.pop("msg", None)
        try:
            if opts.alpha == "custom":
                opts.custom_alpha = d.number(text("customAlpha"), "Alpha")
                if opts.custom_alpha is None or opts.custom_alpha <= 0:
                    raise ValueError("Enter a positive custom alpha.")
            if text("centre").strip():
                opts.centre_points = _whole(text("centre"), "Centre points", 0)
            if text("axialCentre").strip():
                opts.axial_centre_points = _whole(text("axialCentre"),
                                                  "Centre points in the axial block", 0)
            opts.replicates = _whole(text("replicates") or "1", "Replicates", 1)
            opts.blocks = _whole(text("blocks") or "1", "Blocks", 1)
            if text("seed").strip():
                opts.seed = _whole(text("seed"), "The random seed", 0)
        except ValueError as exc:
            problems["msg"] = str(exc)
        return opts

    return d.run(ctx, doc, parent, dialog, title, read, lambda _: problems.get("msg"))


MODELS = (("2fi", "Main effects and two-factor interactions"), ("main", "Main effects only"))


def show_screening_analysis(ctx, doc, parent):
    """Show the dialog; returns ScreeningOptions, or None if cancelled."""
    title = "LibreStats - Analyse Screening Design"
    detected = default_table_range(doc)
    range_text = detected or d.initial_range(doc)[0]
    names, _, problem = response_names(doc, range_text)
    if detected and problem:
        # the active sheet is a Design sheet that cannot be analysed: say why, at once
        calc.message_box(ctx, parent, title, problem)
        return None
    responses = tuple((n, n) for n in names) or \
        (("", "(select a Design table to list its responses)"),)
    b = d.Builder(ctx, title, 260, 164)
    b.label("lblRange", 10, 8, 55, "Design table:")
    b.edit("range", 70, 8, 180, range_text)
    b.label("lblHint", 10, 22, 240, "A Design sheet made with New Design, responses filled in.")
    b.label("lblResponse", 10, 40, 55, "Response:")
    b.dropdown("response", 70, 40, 180, responses)

    b.line("lineModel", 58, "Model")
    b.radio("model2fi", 10, 70, 240, MODELS[0][1], True)
    b.radio("modelMain", 10, 84, 240, MODELS[1][1], False)
    b.label("lblAlpha", 10, 102, 85, "Significance level:")
    b.edit("alpha", 100, 102, 40, "0.05")
    b.check("charts", 10, 120, 240, "Pareto chart and half-normal plot", True)
    b.ok_cancel(144)

    dialog = b.create(parent)
    d.select(dialog, "response", responses, responses[0][0])
    problems = {}

    def read(dlg):
        text = lambda name: dlg.getControl(name).getText()  # noqa: E731
        opts = ScreeningOptions(
            range_text=text("range"), response=d.picked(dlg, "response", responses),
            model="2fi" if d.checked(dlg, "model2fi") else "main",
            charts=d.checked(dlg, "charts"))
        problems.pop("msg", None)
        try:
            opts.alpha = d.number(text("alpha"), "The significance level")
            if opts.alpha is None or not 0 < opts.alpha < 1:
                raise ValueError("The significance level must be between 0 and 1.")
        except ValueError as exc:
            problems["msg"] = str(exc)
        return opts

    return d.run(ctx, doc, parent, dialog, title, read, lambda _: problems.get("msg"))


RSM_MODELS = (("quadratic", "modelQuadratic", "Full quadratic (linear, squares, interactions)"),
              ("interaction", "modelInteraction", "Linear and two-factor interactions"),
              ("linear", "modelLinear", "Linear"))


def show_rsm_analysis(ctx, doc, parent):
    """Show the dialog; returns RSMOptions, or None if cancelled."""
    title = "LibreStats - Analyse Response Surface Design"
    detected = default_table_range(doc)
    range_text = detected or d.initial_range(doc)[0]
    names, factors, problem = response_names(doc, range_text, screening=False)
    if detected and problem:
        calc.message_box(ctx, parent, title, problem)
        return None
    responses = tuple((n, n) for n in names) or \
        (("", "(select a Design table to list its responses)"),)
    factor_items = tuple((f, f) for f in factors) or (("", "(first two factors)"),)
    b = d.Builder(ctx, title, 260, 226)
    b.label("lblRange", 10, 8, 55, "Design table:")
    b.edit("range", 70, 8, 180, range_text)
    b.label("lblHint", 10, 22, 240, "A Design sheet made with New Design, responses filled in.")
    b.label("lblResponse", 10, 40, 55, "Response:")
    b.dropdown("response", 70, 40, 180, responses)

    b.line("lineModel", 58, "Model")
    for i, (_, name, label) in enumerate(RSM_MODELS):
        b.radio(name, 10, 70 + 14 * i, 240, label, i == 0)
    b.label("lblAlpha", 10, 116, 85, "Significance level:")
    b.edit("alpha", 100, 116, 40, "0.05")

    b.line("lineMap", 134, "Response map and 3D surface (other factors at their centre)")
    b.check("charts", 10, 146, 240, "Draw the response map and the 3D surface", True)
    b.label("lblX", 10, 162, 85, "Horizontal factor:")
    b.dropdown("xFactor", 100, 162, 150, factor_items)
    b.label("lblY", 10, 178, 85, "Vertical factor:")
    b.dropdown("yFactor", 100, 178, 150, factor_items)
    b.ok_cancel(206)

    dialog = b.create(parent)
    d.select(dialog, "response", responses, responses[0][0])
    d.select(dialog, "xFactor", factor_items, factor_items[0][0])
    d.select(dialog, "yFactor", factor_items, factor_items[min(1, len(factor_items) - 1)][0])
    problems = {}

    def read(dlg):
        text = lambda name: dlg.getControl(name).getText()  # noqa: E731
        opts = RSMOptions(
            range_text=text("range"), response=d.picked(dlg, "response", responses),
            model=next(m for m, name, _ in RSM_MODELS if d.checked(dlg, name)),
            charts=d.checked(dlg, "charts"),
            x_factor=d.picked(dlg, "xFactor", factor_items),
            y_factor=d.picked(dlg, "yFactor", factor_items))
        problems.pop("msg", None)
        try:
            opts.alpha = d.number(text("alpha"), "The significance level")
            if opts.alpha is None or not 0 < opts.alpha < 1:
                raise ValueError("The significance level must be between 0 and 1.")
            if opts.charts and opts.x_factor and opts.x_factor == opts.y_factor:
                raise ValueError("Choose two different factors for the response map.")
        except ValueError as exc:
            problems["msg"] = str(exc)
        return opts

    return d.run(ctx, doc, parent, dialog, title, read, lambda _: problems.get("msg"))


DIALOGS = {
    "doeAnalyseRSM": show_rsm_analysis,
    "doeAnalyseScreening": show_screening_analysis,
    "doeCCD": lambda ctx, doc, parent: show_rsm(ctx, doc, parent, "ccd"),
    "doeBoxBehnken": lambda ctx, doc, parent: show_rsm(ctx, doc, parent, "bbd"),
    "doeFull": lambda ctx, doc, parent: show_design(ctx, doc, parent, "full"),
    "doeFractional": lambda ctx, doc, parent: show_design(ctx, doc, parent, "fractional"),
    "doePlackettBurman": lambda ctx, doc, parent: show_design(ctx, doc, parent, "pb"),
}
