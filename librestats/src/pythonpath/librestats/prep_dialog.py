"""Dialogs for transformations, scaling, log-ratios and (re)coding."""

import unohelper
from com.sun.star.awt import XActionListener

from . import calc, coding
from . import dialogs as d
from .prep_calc import PrepOptions

BASES = (("e", "e (natural log, ln)"), ("10", "10"), ("2", "2"))
OFFSETS = (("none", "none"), ("one", "+1"), ("half_min", "+ half the smallest positive value"),
           ("custom", "custom:"))
K_BINS = tuple((k, str(k)) for k in range(2, 11))
TIES = (("average", "average (as RANK.AVG)"), ("min", "minimum (as RANK.EQ)"),
        ("max", "maximum"), ("dense", "dense (1, 2, 2, 3)"))
SCORES = (("blom", "Blom (r - 3/8) / (n + 1/4)"), ("vdw", "van der Waerden r / (n + 1)"),
          ("tukey", "Tukey (r - 1/3) / (n + 1/3)"), ("rankit", "Rankit (r - 1/2) / n"))


def _start(ctx, doc, title, height):
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, height)
    b.range_and_labels(range_text, labels)
    return b


def _common(dialog, settings):
    return PrepOptions(dialog.getControl("range").getText(), d.checked(dialog, "labels"),
                       d.checked(dialog, "outBeside"), settings)


def _run(ctx, doc, parent, b, title, read):
    """Create and execute; a ValueError from read() is shown and the dialog reopened."""
    dialog = b.create(parent)
    problems = {}

    def safe_read(dlg):
        try:
            problems.pop("msg", None)
            return read(dlg)
        except ValueError as exc:
            problems["msg"] = str(exc)
            return PrepOptions(dlg.getControl("range").getText())

    return dialog, (lambda: d.run(ctx, doc, parent, dialog, title, safe_read,
                                  lambda _: problems.get("msg")))


# --- Transform -----------------------------------------------------------------------

def show_transform(ctx, doc, parent):
    title = "LibreStats - Transform"
    b = _start(ctx, doc, title, 248)
    b.line("lineMethod", 42, "Transformation")
    methods = ("log", "sqrt", "boxcox", "yeojohnson", "logit")
    b.radio("log", 10, 54, 120, "Logarithm, base:", True)
    b.radio("sqrt", 10, 68, 240, "Square root", False)
    b.radio("boxcox", 10, 82, 240, "Box-Cox (positive values)", False)
    b.radio("yeojohnson", 10, 96, 240, "Yeo-Johnson (any values, also negative)", False)
    b.radio("logit", 10, 110, 240, "Logit  ln(p / (1 - p))  for proportions", False)
    b.dropdown("base", 135, 54, 115, BASES)
    b.label("lblOffset", 10, 128, 120, "Offset (log, square root, Box-Cox):")
    b.dropdown("offset", 135, 128, 75, OFFSETS)
    b.edit("offsetValue", 214, 128, 36, "")
    b.label("lblLambda", 10, 144, 120, "Lambda (Box-Cox, Yeo-Johnson):")
    b.edit("lambda", 135, 144, 40, "")
    b.label("lblLambdaHint", 180, 144, 70, "empty = estimate")
    b.check("percent", 10, 160, 240, "Logit: values are percentages (0-100)", False)
    b.check("squeeze", 10, 172, 240, "Logit: squeeze values of exactly 0 or 1 "
            "(Smithson & Verkuilen, 2006)", True)
    b.output_choice(188)
    b.ok_cancel(230)

    def read(dlg):
        offset = d.picked(dlg, "offset", OFFSETS)
        if offset == "custom":
            offset = d.number(dlg.getControl("offsetValue").getText(), "The custom offset")
            if offset is None:
                raise ValueError("Enter a custom offset value.")
        lam = d.number(dlg.getControl("lambda").getText(), "Lambda")
        return _common(dlg, {"method": d.radio_choice(dlg, methods),
                             "base": d.picked(dlg, "base", BASES), "offset": offset,
                             "lambda": "" if lam is None else lam,
                             "percent": d.checked(dlg, "percent"),
                             "squeeze": d.checked(dlg, "squeeze")})

    dialog, go = _run(ctx, doc, parent, b, title, read)
    d.select(dialog, "base", BASES, "e")
    d.select(dialog, "offset", OFFSETS, "none")
    return go()


# --- Scale ---------------------------------------------------------------------------

def show_scale(ctx, doc, parent):
    title = "LibreStats - Scale"
    b = _start(ctx, doc, title, 160)
    b.line("lineMethod", 42, "Scaling (van den Berg et al., 2006)")
    methods = ("centre", "z", "range", "pareto")
    b.radio("centre", 10, 54, 240, "Centring only:  x - mean", False)
    b.radio("z", 10, 68, 240, "Autoscaling (z-score):  (x - mean) / SD", True)
    b.radio("range", 10, 82, 240, "Range (0-1):  (x - min) / (max - min)", False)
    b.radio("pareto", 10, 96, 240, "Pareto:  (x - mean) / sqrt(SD)", False)
    b.output_choice(112)
    b.ok_cancel(142)
    return _run(ctx, doc, parent, b, title,
                lambda dlg: _common(dlg, {"method": d.radio_choice(dlg, methods)}))[1]()


# --- Compositional -------------------------------------------------------------------

def show_logratio(ctx, doc, parent):
    title = "LibreStats - Compositional Log-ratios"
    b = _start(ctx, doc, title, 176)
    b.line("lineMethod", 42, "Log-ratio transformation (Aitchison, 1982)")
    b.label("lblInfo", 10, 52, 240, "Each row is one composition; every selected column is a part.")
    b.radio("clr", 10, 66, 240, "Centred log-ratio (CLR): one column per part", True)
    b.radio("ilr", 10, 80, 240, "Isometric log-ratio (ILR, pivot coordinates): D - 1 columns",
            False)
    b.check("replace", 10, 98, 240, "Replace zeros (multiplicative, 0.65 x smallest positive "
            "value)", False)
    b.output_choice(114)
    b.ok_cancel(156)
    return _run(ctx, doc, parent, b, title, lambda dlg: _common(dlg, {
        "method": d.radio_choice(dlg, ("clr", "ilr")),
        "zeros": "replace" if d.checked(dlg, "replace") else "error"}))[1]()


# --- Coding --------------------------------------------------------------------------

def show_dummies(ctx, doc, parent):
    title = "LibreStats - Dummy Variables"
    b = _start(ctx, doc, title, 134)
    b.line("lineMethod", 42, "One 0/1 column per category of each selected column")
    b.check("dropFirst", 10, 56, 240, "Leave out the first category (reference; k - 1 columns "
            "for regression)", False)
    b.output_choice(74)
    b.ok_cancel(114)
    return _run(ctx, doc, parent, b, title,
                lambda dlg: _common(dlg, {"drop_first": d.checked(dlg, "dropFirst")}))[1]()


def show_presence(ctx, doc, parent):
    title = "LibreStats - Presence/Absence"
    b = _start(ctx, doc, title, 146)
    b.line("lineMethod", 42, "1 = present, 0 = absent")
    b.label("lblThreshold", 10, 56, 120, "Present if the value is above:")
    b.edit("threshold", 135, 56, 40, "0")
    b.check("emptyAbsent", 10, 72, 240, "Empty cells mean absent (0) instead of missing", False)
    b.output_choice(88)
    b.ok_cancel(128)

    def read(dlg):
        threshold = d.number(dlg.getControl("threshold").getText(), "The threshold")
        return _common(dlg, {"threshold": 0.0 if threshold is None else threshold,
                             "empty_absent": d.checked(dlg, "emptyAbsent")})

    return _run(ctx, doc, parent, b, title, read)[1]()


def show_bins(ctx, doc, parent):
    title = "LibreStats - Discretize into Bins"
    b = _start(ctx, doc, title, 160)
    b.line("lineMethod", 42, "Bin boundaries")
    methods = ("width", "quantile", "kmeans")
    b.radio("width", 10, 54, 240, "Equal width (between minimum and maximum)", True)
    b.radio("quantile", 10, 68, 240, "Quantiles (about equal counts per bin)", False)
    b.radio("kmeans", 10, 82, 240, "k-means (optimal grouping of similar values)", False)
    b.label("lblK", 10, 98, 60, "Number of bins:")
    b.dropdown("k", 75, 98, 40, K_BINS)
    b.output_choice(114)
    b.ok_cancel(142)
    dialog, go = _run(ctx, doc, parent, b, title, lambda dlg: _common(dlg, {
        "method": d.radio_choice(dlg, methods), "k": d.picked(dlg, "k", K_BINS)}))
    d.select(dialog, "k", K_BINS, 4)
    return go()


class _LoadCategories(unohelper.Base, XActionListener):
    """'List categories' button: fill the rules box with 'value -> value' lines."""

    def __init__(self, ctx, doc, parent, dialog):
        self.ctx, self.doc, self.parent, self.dialog = ctx, doc, parent, dialog

    def actionPerformed(self, event):
        try:
            rng = calc.resolve_range(self.doc, self.dialog.getControl("range").getText())
            data = calc.used_data(rng)
            if d.checked(self.dialog, "labels"):
                data = data[1:]
            values = [v for row in data for v in row]
            self.dialog.getControl("rules").setText(coding.rules_template(values))
        except calc.RangeError as exc:
            calc.message_box(self.ctx, self.parent, "LibreStats", str(exc))

    def disposing(self, event):
        pass


def show_recode(ctx, doc, parent):
    title = "LibreStats - Recode Categories"
    b = _start(ctx, doc, title, 236)
    b.line("lineRules", 42, "Rules: old value(s) -> new value")
    b.label("lblSyntax", 10, 52, 240, "Separate several old values with ';'.  * -> x  recodes "
            "everything else.")
    b.multiline("rules", 10, 66, 240, 100, "")
    b.button("load", 10, 170, 110, "List categories in the range")
    b.label("lblKeep", 125, 172, 125, "Unlisted values are kept.")
    b.output_choice(188)
    b.ok_cancel(218)
    dialog, go = _run(ctx, doc, parent, b, title, lambda dlg: _common(
        dlg, {"rules": dlg.getControl("rules").getText()}))
    listener = _LoadCategories(ctx, doc, parent, dialog)
    dialog.getControl("load").addActionListener(listener)
    listener.actionPerformed(None)  # start with the categories of the current selection
    return go()


def show_ranks(ctx, doc, parent):
    title = "LibreStats - Ranks and Normal Scores"
    b = _start(ctx, doc, title, 150)
    b.line("lineMethod", 42, "New columns")
    b.check("ranks", 10, 56, 90, "Ranks, ties:", True)
    b.dropdown("ties", 105, 56, 145, TIES)
    b.check("scores", 10, 72, 90, "Normal scores:", False)
    b.dropdown("scoreMethod", 105, 72, 145, SCORES)
    b.output_choice(90)
    b.ok_cancel(130)

    def read(dlg):
        return _common(dlg, {"ranks": d.checked(dlg, "ranks"), "ties": d.picked(dlg, "ties", TIES),
                             "scores": d.checked(dlg, "scores"),
                             "score_method": d.picked(dlg, "scoreMethod", SCORES)})

    dialog, go = _run(ctx, doc, parent, b, title, read)
    d.select(dialog, "ties", TIES, "average")
    d.select(dialog, "scoreMethod", SCORES, "blom")
    return go()


DIALOGS = {"transform": show_transform, "scale": show_scale, "logratio": show_logratio,
           "dummies": show_dummies, "presence": show_presence, "bins": show_bins,
           "recode": show_recode, "ranks": show_ranks}
