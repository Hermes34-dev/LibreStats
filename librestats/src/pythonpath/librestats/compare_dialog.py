"""Dialogs for Comparison Tests > t-Tests / ANOVA."""

from . import calc
from . import dialogs as d
from .assumptions_dialog import ALPHAS, _groups_and_alpha
from .compare_calc import AnovaOptions, TTestOptions

ALTERNATIVES = (("two-sided", "Two-sided (≠)"), ("greater", "One-sided, greater (>)"),
                ("less", "One-sided, less (<)"))


def _run(ctx, doc, parent, b, title, make, groups=True, defaults=()):
    """Show the dialog until the input is valid; make(dlg) builds the options.
    defaults: further (dropdown name, items, value) to preselect."""
    dialog = b.create(parent)
    d.select(dialog, "alpha", ALPHAS, 0.05)
    if dialog.getControl("alternative"):
        d.select(dialog, "alternative", ALTERNATIVES, "two-sided")
    for name, items, value in defaults:
        d.select(dialog, name, items, value)
    problems = {}

    def read(dlg):
        problems.pop("msg", None)
        group_text = dlg.getControl("groupRange").getText() if groups else ""
        try:
            opts = make(dlg)
        except ValueError as exc:
            problems["msg"] = str(exc)
            return TTestOptions(dlg.getControl("range").getText())
        opts.range_text = dlg.getControl("range").getText()
        opts.has_labels = d.checked(dlg, "labels")
        opts.alpha = d.picked(dlg, "alpha", ALPHAS)
        if groups:
            opts.group_range = group_text
            if group_text.strip():
                try:
                    calc.resolve_range(doc, group_text)
                except calc.RangeError as exc:
                    problems["msg"] = "Group labels: " + str(exc)
        return opts

    return d.run(ctx, doc, parent, dialog, title, read, lambda _: problems.get("msg"))


def _alternative(dlg):
    return d.picked(dlg, "alternative", ALTERNATIVES)


def show_one_sample(ctx, doc, parent):
    title = "LibreStats - One-Sample t-Test"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 172)
    b.range_and_labels(range_text, labels)
    _groups_and_alpha(b, 42, "Empty: test each column. Else: one column of values, per group.")
    b.line("lineTest", 104, "Test")
    b.label("lblMu", 10, 116, 85, "Test value (μ0):")
    b.edit("mu", 100, 116, 60, "0")
    b.label("lblAlt", 10, 132, 85, "Alternative:")
    b.dropdown("alternative", 100, 132, 150, ALTERNATIVES)
    b.ok_cancel(152)

    def make(dlg):
        mu = d.number(dlg.getControl("mu").getText(), "The test value")
        return TTestOptions(mu=0.0 if mu is None else mu, alternative=_alternative(dlg))

    return _run(ctx, doc, parent, b, title, make)


def show_welch(ctx, doc, parent):
    title = "LibreStats - Welch's t-Test"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 172)
    b.range_and_labels(range_text, labels)
    _groups_and_alpha(b, 42, "Empty: two columns. Else: one column of values, two groups.")
    b.line("lineTest", 104, "Test (difference = first - second sample)")
    b.label("lblAlt", 10, 116, 85, "Alternative:")
    b.dropdown("alternative", 100, 116, 150, ALTERNATIVES)
    b.label("lblInfo", 10, 132, 240, "Does not assume equal variances (Welch, 1947).")
    b.ok_cancel(152)
    return _run(ctx, doc, parent, b, title,
                lambda dlg: TTestOptions(alternative=_alternative(dlg)))


def show_paired(ctx, doc, parent):
    title = "LibreStats - Paired t-Test"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 130)
    b.range_and_labels(range_text, labels)
    b.line("lineTest", 42, "Test (two columns, one row per pair; difference = first - second)")
    b.label("lblAlt", 10, 54, 85, "Alternative:")
    b.dropdown("alternative", 100, 54, 150, ALTERNATIVES)
    b.label("lblAlpha", 10, 70, 85, "Significance level:")
    b.dropdown("alpha", 100, 70, 50, ALPHAS)
    b.label("lblInfo", 10, 88, 240, "Rows with a missing value are left out.")
    b.ok_cancel(110)
    return _run(ctx, doc, parent, b, title,
                lambda dlg: TTestOptions(alternative=_alternative(dlg)), groups=False)


def _show_anova(ctx, doc, parent, title, info):
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 150)
    b.range_and_labels(range_text, labels)
    _groups_and_alpha(b, 42, "Empty: each column is a group. Else: one column of values.")
    b.label("lblInfo", 10, 108, 240, info)
    b.ok_cancel(130)
    return _run(ctx, doc, parent, b, title, lambda dlg: AnovaOptions())


def show_anova(ctx, doc, parent):
    return _show_anova(ctx, doc, parent, "LibreStats - One-Way ANOVA",
                       "Classical F-test; assumes equal variances in all groups.")


def show_welch_anova(ctx, doc, parent):
    return _show_anova(ctx, doc, parent, "LibreStats - Welch's ANOVA",
                       "Welch's F-test (1951); does not assume equal variances.")


def show_rm_anova(ctx, doc, parent):
    title = "LibreStats - Repeated-Measures ANOVA"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 124)
    b.range_and_labels(range_text, labels)
    b.line("lineTest", 42, "Test (one column per condition, one row per subject)")
    b.label("lblAlpha", 10, 54, 85, "Significance level:")
    b.dropdown("alpha", 100, 54, 50, ALPHAS)
    b.label("lblInfo", 10, 72, 240, "Rows with a missing value are left out. Reports Mauchly's")
    b.label("lblInfo2", 10, 84, 240, "test and Greenhouse-Geisser / Huynh-Feldt corrections.")
    b.ok_cancel(104)
    return _run(ctx, doc, parent, b, title, lambda dlg: AnovaOptions(), groups=False)


DIALOGS = {"ttestOne": show_one_sample, "ttestWelch": show_welch, "ttestPaired": show_paired,
           "anovaOneWay": show_anova, "anovaWelch": show_welch_anova, "anovaRM": show_rm_anova}
