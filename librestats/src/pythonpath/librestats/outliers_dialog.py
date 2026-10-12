"""Options dialog for Data Preparation > Outliers."""

from . import dialogs as d
from .outliers_calc import OutlierOptions

TUKEY_K = ((1.5, "1.5 (Tukey)"), (2.2, "2.2 (Hoaglin & Iglewicz)"), (3.0, "3.0 (far out)"))
MAD_THRESHOLDS = ((2.0, "2.0 (lenient)"), (2.5, "2.5 (recommended)"), (3.0, "3.0 (strict)"))
ALPHAS = ((0.10, "0.10"), (0.05, "0.05"), (0.01, "0.01"))
MAHALANOBIS_ALPHAS = ((0.001, "0.001 (usual)"), (0.01, "0.01"), (0.05, "0.05"))

TITLE = "LibreStats - Outliers"


def show(ctx, doc, parent):
    """Show the dialog; returns OutlierOptions, or None if cancelled."""
    opts = OutlierOptions()
    opts.range_text, opts.has_labels = d.initial_range(doc)

    b = d.Builder(ctx, TITLE, 260, 236)
    b.range_and_labels(opts.range_text, opts.has_labels)

    b.line("lineUni", 42, "Univariate (each column tested separately)")
    b.check("tukey", 10, 56, 135, "Tukey fences (IQR)", opts.tukey)
    b.label("lblK", 150, 56, 15, "k:")
    b.dropdown("tukeyK", 165, 56, 85, TUKEY_K)
    b.check("mad", 10, 72, 135, "Robust z-scores (median / MAD)", opts.mad)
    b.label("lblZ", 150, 72, 15, "|z| >")
    b.dropdown("madThreshold", 165, 72, 85, MAD_THRESHOLDS)
    b.check("grubbs", 10, 88, 135, "Grubbs test", opts.grubbs)
    b.check("grubbsIterate", 150, 88, 100, "Repeat until none found", opts.grubbs_iterate)
    b.check("dixon", 10, 104, 135, "Dixon test (3 to 30 values)", opts.dixon)
    b.label("lblAlpha", 10, 122, 135, "Significance level (Grubbs, Dixon):")
    b.dropdown("alpha", 165, 122, 85, ALPHAS)

    b.line("lineMulti", 142, "Multivariate (all columns together)")
    b.check("mahalanobis", 10, 156, 135, "Mahalanobis distance", opts.mahalanobis)
    b.label("lblMahalAlpha", 150, 156, 15, "alpha:")
    b.dropdown("mahalanobisAlpha", 165, 156, 85, MAHALANOBIS_ALPHAS)

    b.line("lineOut", 176, "Output")
    b.label("lblOut", 10, 188, 240, "Results are written to a new sheet.")
    b.check("highlight", 10, 200, 240, "Highlight outliers in the input data", opts.highlight)
    b.ok_cancel(218)

    dialog = b.create(parent)
    d.select(dialog, "tukeyK", TUKEY_K, opts.tukey_k)
    d.select(dialog, "madThreshold", MAD_THRESHOLDS, opts.mad_threshold)
    d.select(dialog, "alpha", ALPHAS, opts.alpha)
    d.select(dialog, "mahalanobisAlpha", MAHALANOBIS_ALPHAS, opts.mahalanobis_alpha)

    def validate(result):
        return None if result.any_method() else "Please select at least one method."

    return d.run(ctx, doc, parent, dialog, TITLE, _read, validate)


def _read(dialog):
    return OutlierOptions(
        range_text=dialog.getControl("range").getText(),
        has_labels=d.checked(dialog, "labels"),
        tukey=d.checked(dialog, "tukey"), tukey_k=d.picked(dialog, "tukeyK", TUKEY_K),
        mad=d.checked(dialog, "mad"),
        mad_threshold=d.picked(dialog, "madThreshold", MAD_THRESHOLDS),
        grubbs=d.checked(dialog, "grubbs"), grubbs_iterate=d.checked(dialog, "grubbsIterate"),
        dixon=d.checked(dialog, "dixon"), alpha=d.picked(dialog, "alpha", ALPHAS),
        mahalanobis=d.checked(dialog, "mahalanobis"),
        mahalanobis_alpha=d.picked(dialog, "mahalanobisAlpha", MAHALANOBIS_ALPHAS),
        highlight=d.checked(dialog, "highlight"))
