"""Dialogs for new screening and response surface designs (DOE / RSM menu)."""

import unohelper
from com.sun.star.awt import XItemListener

from . import calc
from . import dialogs as d
from . import doe
from . import rsm
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


RSM_KIND_RADIOS = (("ccd", "kindCCD"), ("bbd", "kindBB"), ("cccd", "kindCCCD"))


def _factor_count(doc, opts):
    """Number of factors in the factor table, or None if it cannot be read yet."""
    try:
        data = calc.used_data(calc.resolve_range(doc, opts.range_text))
        return len(doe.read_factors(data[1:] if opts.has_labels else data))
    except Exception:                    # the range check and run_design report problems
        return None


def show_rsm(ctx, doc, parent, kind="ccd"):
    """Show the dialog; returns DesignOptions, or None if cancelled."""
    title = "LibreStats - New Response Surface Design"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 390)
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
    b.radio("kindCCCD", 10, 94, 240, "CCCD (Prieto et al.): CCD + outer cube at +-alpha, "
            "5 levels", kind == "cccd")
    b.label("lblCCCD1", 20, 106, 230, "CCCD: Low/High = extreme levels (+-alpha); every run "
            "stays within them.")
    b.label("lblCCCD2", 20, 116, 230, "Leichtweis, Pereira, Prieto et al. (2019) Molecules "
            "24(3):573. Most textbooks")
    b.label("lblCCCD3", 20, 126, 230, "call the plain 20-run CCD 'circumscribed'; this "
            "variant adds the outer cube.")

    y = 36   # offset of everything below the CCCD notes
    b.line("lineCCD", 96 + y, "Axial distance (CCCD: (2^k)^(1/4) or custom)")
    b.label("lblAlpha", 10, 108 + y, 85, "Axial distance alpha:")
    b.radio("alphaRotatable", 100, 108 + y, 150, "rotatable (CCCD: (2^k)^(1/4))", True)
    b.radio("alphaOrthogonal", 100, 122 + y, 150, "orthogonal blocking (needs blocks)", False)
    b.radio("alphaFace", 100, 136 + y, 150, "face-centred (alpha = 1)", False)
    b.radio("alphaCustom", 100, 150 + y, 50, "custom:", False)
    b.edit("customAlpha", 150, 150 + y, 40, "")
    b.check("fractionalCube", 10, 166 + y, 240, "CCD: resolution V fraction as factorial "
            "part (5 or more factors)", True)
    b.check("inscribed", 10, 180 + y, 240, "CCD inscribed: axial points at Low and High "
            "(limits that cannot be exceeded)", False)

    b.line("lineOptions", 196 + y, "Options (leave centre points empty for the usual numbers)")
    b.label("lblCentre", 10, 208 + y, 120, "Centre points (per cube block):")
    b.edit("centre", 135, 208 + y, 40, "")
    b.label("lblAxialCentre", 10, 224 + y, 120, "Centre points in the axial block:")
    b.edit("axialCentre", 135, 224 + y, 40, "")
    b.label("lblReplicates", 10, 240 + y, 120, "Replicates:")
    b.edit("replicates", 135, 240 + y, 40, "1")
    b.label("lblBlocks", 10, 256 + y, 120, "Blocks:")
    b.edit("blocks", 135, 256 + y, 40, "1")
    b.label("lblBlocksHint", 10, 270 + y, 240, "Replicates x 1, 2, 3 or 5 (CCD); x 1, or 3 "
            "for a 4-factor Box-Behnken; x 1 (CCCD).")
    b.check("randomize", 10, 286 + y, 120, "Randomise the run order", True)
    b.label("lblSeed", 130, 286 + y, 40, "Seed:")
    b.edit("seed", 170, 286 + y, 50, "")

    b.line("lineResponses", 302 + y, "Responses")
    b.label("lblResponses", 10, 314 + y, 85, "Response columns:")
    b.edit("responses", 100, 314 + y, 150, "Response")
    b.ok_cancel(334 + y)

    dialog = b.create(parent)
    problems = {}

    def read(dlg):
        text = lambda name: dlg.getControl(name).getText()  # noqa: E731
        opts = DesignOptions(
            range_text=text("range"), has_labels=d.checked(dlg, "labels"),
            kind=next(k for k, name in RSM_KIND_RADIOS if d.checked(dlg, name)),
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
            if opts.kind == "cccd" and opts.alpha in ("orthogonal", "face"):
                raise ValueError("The CCCD uses alpha = (2^k)^(1/4) (the first choice) or a "
                                 "custom alpha greater than 1.")
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

    def validate(opts):
        if problems.get("msg"):
            return problems["msg"]
        if opts.kind == "cccd":
            k = _factor_count(doc, opts)
            if k is not None and k >= doe.CCCD_LARGE:
                runs = (2 * 2 ** k + 2 * k) * opts.replicates
                if not calc.confirm(ctx, parent, title,
                                    f"The outer cube doubles the factorial part: with {k} "
                                    f"factors the CCCD has {runs} runs plus centre points. "
                                    "Generate it anyway?"):
                    return ""                # back to the dialog
        return None

    return d.run(ctx, doc, parent, dialog, title, read, validate)


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
              ("linear", "modelLinear", "Linear"),
              ("cubic", "modelCubic", "Reduced cubic (full cubic, aliased terms removed)"))


REDUCTION_RADIOS = (("none", "redNone", "None: fit the chosen model as is"),
                    ("manual", "redManual", "Manual: remove the terms selected below"),
                    ("backward", "redBackward", "Backward elimination (p-value), alpha-out:"),
                    ("subsets", "redSubsets", "All subsets (AICc): show top"))
ALPHA_OUT = (("0.05", "0.05"), ("0.10", "0.10"))
FIXED_RADIOS = (("optimum", "fixOpt", "optimum", 42), ("centre", "fixCentre", "centre", 37),
                ("custom", "fixCustom", "custom:", 38))


class _TermLists(unohelper.Base, XItemListener):
    """Fill the 'remove' and 'locked' lists with the terms of the chosen model."""

    def __init__(self, dialog, k):
        self.dialog, self.k = dialog, k

    def update(self):
        model = next((m for m, name, _ in RSM_MODELS if d.checked(self.dialog, name)),
                     "quadratic")
        names = tuple(t.name for t in rsm.model_terms(self.k, model)[1:]) if self.k else ()
        for name in ("removeTerms", "lockedTerms"):
            self.dialog.getControl(name).getModel().StringItemList = names

    def itemStateChanged(self, event):
        self.update()

    def disposing(self, event):
        pass


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
    b = d.Builder(ctx, title, 260, 426)
    b.label("lblRange", 10, 8, 55, "Design table:")
    b.edit("range", 70, 8, 180, range_text)
    b.label("lblHint", 10, 22, 240, "A Design sheet made with New Design, responses filled in.")
    b.label("lblResponse", 10, 40, 55, "Response:")
    b.dropdown("response", 70, 40, 180, responses)

    b.line("lineModel", 58, "Model")
    for i, (_, name, label) in enumerate(RSM_MODELS):
        b.radio(name, 10, 70 + 14 * i, 240, label, i == 0)
    b.label("lblAlpha", 10, 130, 85, "Significance level:")
    b.edit("alpha", 100, 130, 40, "0.05")

    b.line("lineReduction", 148, "Model reduction")
    for i, (_, name, label) in enumerate(REDUCTION_RADIOS):
        b.radio(name, 10, 160 + 14 * i, 165, label, i == 0)
    b.dropdown("alphaOut", 178, 188, 40, ALPHA_OUT)
    b.edit("top", 178, 202, 22, "5")
    b.label("lblRank", 203, 202, 25, "use rank")
    b.edit("rank", 230, 202, 20, "1")
    b.check("hierarchy", 10, 218, 240, "Keep model hierarchy (a term keeps its lower-order "
            "terms)", True)
    b.check("excludeAliased", 10, 232, 240, "Reduced cubic: exclude aliased terms from "
            "selection (e.g. AB^2 = AC^2)", False)
    b.label("lblRemove", 10, 246, 115, "Terms to remove (Manual):")
    b.label("lblLocked", 135, 246, 115, "Locked terms (never removed):")
    b.listbox("removeTerms", 10, 260, 115, 52)
    b.listbox("lockedTerms", 135, 260, 115, 52)
    b.label("lblTermsHint", 10, 314, 240, "Ctrl+click selects several. The constant and "
            "block terms are always kept.")

    b.line("lineMap", 330, "Charts (maps and surfaces: every factor pair, at most 3)")
    b.check("charts", 10, 342, 240, "Maps, 3D surfaces, observed vs predicted, residuals by "
            "factor, optimum trace, influence plot", True)
    b.label("lblFix", 10, 356, 70, "Other factors fixed at:")
    for i, (_, name, label, w) in enumerate(FIXED_RADIOS):
        b.radio(name, (82, 125, 163)[i], 356, w, label, i == 0)
    b.edit("fixedValues", 203, 356, 47, "")
    b.label("lblGoal", 10, 370, 38, "Optimum:")
    b.radio("goalMax", 50, 370, 42, "maximum", True)
    b.radio("goalMin", 93, 370, 42, "minimum", False)
    b.label("lblAnchor", 140, 370, 40, "Trace from:")
    b.radio("anchorOpt", 182, 370, 37, "optimum", True)
    b.radio("anchorCentre", 220, 370, 32, "centre", False)
    b.label("lblRes", 10, 384, 38, "Residuals:")
    b.radio("resRaw", 50, 384, 25, "raw", True)
    b.radio("resStud", 77, 384, 70, "ext. studentized", False)
    b.check("greyscale", 160, 384, 90, "Greyscale (journal style)", False)
    b.ok_cancel(404)

    dialog = b.create(parent)
    d.select(dialog, "response", responses, responses[0][0])
    d.select(dialog, "alphaOut", ALPHA_OUT, "0.05")
    dialog.getControl("fixedValues").getModel().HelpText = \
        "Custom: real values, e.g. Time=13; EhOH=50 (names or letters)"
    lists = _TermLists(dialog, len(factors))
    lists.update()
    for _, name, _ in RSM_MODELS:
        dialog.getControl(name).addItemListener(lists)
    problems = {}

    def read(dlg):
        text = lambda name: dlg.getControl(name).getText()  # noqa: E731
        opts = RSMOptions(
            range_text=text("range"), response=d.picked(dlg, "response", responses),
            model=next(m for m, name, _ in RSM_MODELS if d.checked(dlg, name)),
            charts=d.checked(dlg, "charts"),
            reduction=next(m for m, name, _ in REDUCTION_RADIOS if d.checked(dlg, name)),
            reduction_alpha=float(d.picked(dlg, "alphaOut", ALPHA_OUT)),
            hierarchy=d.checked(dlg, "hierarchy"),
            exclude_aliased=d.checked(dlg, "excludeAliased"),
            fixed=next(f for f, name, _, _ in FIXED_RADIOS if d.checked(dlg, name)),
            fixed_values=text("fixedValues"),
            goal="max" if d.checked(dlg, "goalMax") else "min",
            anchor="optimum" if d.checked(dlg, "anchorOpt") else "centre",
            residuals="raw" if d.checked(dlg, "resRaw") else "studentized",
            greyscale=d.checked(dlg, "greyscale"),
            locked=tuple(dlg.getControl("lockedTerms").getSelectedItems()),
            remove=tuple(dlg.getControl("removeTerms").getSelectedItems()))
        problems.pop("msg", None)
        try:
            opts.alpha = d.number(text("alpha"), "The significance level")
            if opts.alpha is None or not 0 < opts.alpha < 1:
                raise ValueError("The significance level must be between 0 and 1.")
            if opts.reduction == "subsets":
                opts.top = _whole(text("top") or "5", "The number of models shown", 1)
                opts.rank = _whole(text("rank") or "1", "The rank of the model used", 1)
            if opts.reduction == "manual" and not opts.remove:
                raise ValueError("Manual model reduction: select the terms to remove in the "
                                 "list 'Terms to remove'.")
            if opts.charts and opts.fixed == "custom" and not opts.fixed_values.strip():
                raise ValueError("Custom fixed values: enter them as name=value, e.g. "
                                 "Time=13; EhOH=50.")
            both = set(opts.remove) & set(opts.locked)
            if opts.reduction == "manual" and both:
                raise ValueError("A term cannot be both removed and locked: "
                                 + ", ".join(sorted(both)) + ".")
        except ValueError as exc:
            problems["msg"] = str(exc)
        return opts

    return d.run(ctx, doc, parent, dialog, title, read, lambda _: problems.get("msg"))


DIALOGS = {
    "doeAnalyseRSM": show_rsm_analysis,
    "doeAnalyseScreening": show_screening_analysis,
    "doeCCD": lambda ctx, doc, parent: show_rsm(ctx, doc, parent, "ccd"),
    "doeBoxBehnken": lambda ctx, doc, parent: show_rsm(ctx, doc, parent, "bbd"),
    "doeCCCD": lambda ctx, doc, parent: show_rsm(ctx, doc, parent, "cccd"),
    "doeFull": lambda ctx, doc, parent: show_design(ctx, doc, parent, "full"),
    "doeFractional": lambda ctx, doc, parent: show_design(ctx, doc, parent, "fractional"),
    "doePlackettBurman": lambda ctx, doc, parent: show_design(ctx, doc, parent, "pb"),
}
