"""Dialogs for Assumption Checks > Normality / Equal Variances."""

from . import calc
from . import dialogs as d
from .assumptions_calc import NormalityOptions, VarianceOptions

ALPHAS = ((0.05, "0.05"), (0.01, "0.01"), (0.10, "0.10"))


def _groups_and_alpha(b, y, hint):
    """Optional group range and significance level (60 units high)."""
    b.line("lineGroups", y, "Groups")
    b.label("lblGroups", 10, y + 12, 85, "Group labels (optional):")
    b.edit("groupRange", 100, y + 12, 150, "")
    b.label("lblGroupsHint", 10, y + 26, 240, hint)
    b.label("lblAlpha", 10, y + 44, 85, "Significance level:")
    b.dropdown("alpha", 100, y + 44, 50, ALPHAS)


def _run(ctx, doc, parent, b, title, make):
    dialog = b.create(parent)
    d.select(dialog, "alpha", ALPHAS, 0.05)
    problems = {}

    def read(dlg):
        groups = dlg.getControl("groupRange").getText()
        opts = make(dlg, dlg.getControl("range").getText(), d.checked(dlg, "labels"), groups,
                    d.picked(dlg, "alpha", ALPHAS))
        problems.pop("msg", None)
        if groups.strip():
            try:
                calc.resolve_range(doc, groups)
            except calc.RangeError as exc:
                problems["msg"] = "Group labels: " + str(exc)
        return opts

    return d.run(ctx, doc, parent, dialog, title, read, lambda _: problems.get("msg"))


def show_normality(ctx, doc, parent):
    title = "LibreStats - Normality (Shapiro-Wilk)"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 150)
    b.range_and_labels(range_text, labels)
    _groups_and_alpha(b, 42, "Empty: test each column. Else: one column of values, per group.")
    b.label("lblInfo", 10, 108, 240, "Shapiro-Wilk test (Royston, 1995), for 3 to 5000 values.")
    b.ok_cancel(130)
    return _run(ctx, doc, parent, b, title,
                lambda dlg, r, lab, g, alpha: NormalityOptions(r, lab, g, alpha))


def show_variances(ctx, doc, parent):
    title = "LibreStats - Equal Variances"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 176)
    b.range_and_labels(range_text, labels)
    _groups_and_alpha(b, 42, "Empty: each column is a group. Else: one column of values.")
    b.line("lineCentre", 104, "Levene's test, deviations from the group")
    b.radio("median", 10, 116, 240, "median: Brown-Forsythe version (recommended)", True)
    b.radio("mean", 10, 130, 240, "mean: original Levene test (assumes normal data)", False)
    b.ok_cancel(156)
    return _run(ctx, doc, parent, b, title,
                lambda dlg, r, lab, g, alpha: VarianceOptions(
                    r, lab, g, "median" if d.checked(dlg, "median") else "mean", alpha))


DIALOGS = {"normality": show_normality, "variances": show_variances}
