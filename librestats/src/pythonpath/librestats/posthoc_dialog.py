"""Dialogs for Comparison Tests > Post Hoc."""

from . import dialogs as d
from .assumptions_dialog import _groups_and_alpha
from .compare_dialog import ALTERNATIVES, _alternative, _run
from .posthoc_calc import PosthocOptions

ADJUSTMENTS = (("holm", "Holm (recommended)"), ("bonferroni", "Bonferroni (Dunn's original)"),
               ("sidak", "Šidák"), ("bh", "Benjamini-Hochberg (false discovery rate)"))

GROUPS_HINT = "Empty: each column is a group. Else: one column of values."


def _show_pairs(ctx, doc, parent, title, info):
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 150)
    b.range_and_labels(range_text, labels)
    _groups_and_alpha(b, 42, GROUPS_HINT)
    b.label("lblInfo", 10, 108, 240, info)
    b.ok_cancel(130)
    return _run(ctx, doc, parent, b, title, lambda dlg: PosthocOptions())


def show_tukey(ctx, doc, parent):
    return _show_pairs(ctx, doc, parent, "LibreStats - Tukey HSD",
                       "All pairs; pooled variance (assumes equal variances).")


def show_games_howell(ctx, doc, parent):
    return _show_pairs(ctx, doc, parent, "LibreStats - Games-Howell",
                       "All pairs; does not assume equal variances.")


def show_dunnett(ctx, doc, parent):
    title = "LibreStats - Dunnett's Test"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 186)
    b.range_and_labels(range_text, labels)
    _groups_and_alpha(b, 42, GROUPS_HINT)
    b.line("lineTest", 104, "Each group against the control (difference = group - control)")
    b.label("lblControl", 10, 116, 85, "Control group:")
    b.edit("control", 100, 116, 150, "")
    b.label("lblControlHint", 10, 130, 240,
            "A column name or group label; empty: the first group.")
    b.label("lblAlt", 10, 146, 85, "Alternative:")
    b.dropdown("alternative", 100, 146, 150, ALTERNATIVES)
    b.ok_cancel(166)
    return _run(ctx, doc, parent, b, title,
                lambda dlg: PosthocOptions(control=dlg.getControl("control").getText(),
                                           alternative=_alternative(dlg)))


def show_dunn(ctx, doc, parent):
    title = "LibreStats - Dunn's Test"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 160)
    b.range_and_labels(range_text, labels)
    _groups_and_alpha(b, 42, GROUPS_HINT)
    b.line("lineTest", 104, "All pairs of mean ranks, after a Kruskal-Wallis test")
    b.label("lblAdjust", 10, 116, 85, "p-value adjustment:")
    b.dropdown("adjust", 100, 116, 150, ADJUSTMENTS)
    b.ok_cancel(140)

    def make(dlg):
        return PosthocOptions(adjust=d.picked(dlg, "adjust", ADJUSTMENTS))

    return _run(ctx, doc, parent, b, title, make, defaults=[("adjust", ADJUSTMENTS, "holm")])


DIALOGS = {"posthocTukey": show_tukey, "posthocGamesHowell": show_games_howell,
           "posthocDunnett": show_dunnett, "posthocDunn": show_dunn}
