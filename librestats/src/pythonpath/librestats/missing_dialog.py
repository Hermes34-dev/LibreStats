"""Dialogs for Data Preparation > Missing Data Diagnostics / Impute Missing Values."""

from . import dialogs as d
from .missing_calc import DiagnoseOptions, ImputeOptions

K_VALUES = ((3, "3"), (5, "5"), (10, "10"))
MICE_METHODS = (("pmm", "Predictive mean matching"), ("norm", "Bayesian linear regression"))
M_VALUES = ((5, "5"), (10, "10"), (20, "20"), (50, "50"))
ITERATIONS = ((5, "5"), (10, "10"), (20, "20"))
METHODS = ("mean", "median", "knn", "mice")

DIAGNOSE_TITLE = "LibreStats - Missing Data Diagnostics"
IMPUTE_TITLE = "LibreStats - Impute Missing Values"


def show_diagnose(ctx, doc, parent):
    opts = DiagnoseOptions()
    opts.range_text, opts.has_labels = d.initial_range(doc)
    b = d.Builder(ctx, DIAGNOSE_TITLE, 260, 104)
    b.range_and_labels(opts.range_text, opts.has_labels)
    b.line("lineTests", 42, "Analyses")
    b.label("lblAlways", 10, 54, 240, "Missing values per variable and missing-data patterns.")
    b.check("little", 10, 68, 240, "Little's MCAR test (are data missing completely at random?)",
            opts.little)
    b.ok_cancel(86)

    def read(dialog):
        return DiagnoseOptions(dialog.getControl("range").getText(),
                               d.checked(dialog, "labels"), d.checked(dialog, "little"))

    return d.run(ctx, doc, parent, b.create(parent), DIAGNOSE_TITLE, read)


def show_impute(ctx, doc, parent):
    opts = ImputeOptions()
    opts.range_text, opts.has_labels = d.initial_range(doc)
    b = d.Builder(ctx, IMPUTE_TITLE, 260, 254)
    b.range_and_labels(opts.range_text, opts.has_labels)

    b.line("lineMethod", 42, "Method")
    # the four method radio buttons are inserted together so they form one group
    b.radio("mean", 10, 56, 240, "Mean of the variable", True)
    b.radio("median", 10, 70, 240, "Median of the variable", False)
    b.radio("knn", 10, 84, 135, "k-nearest neighbours", False)
    b.radio("mice", 10, 116, 240, "Multiple imputation by chained equations (MICE)", False)
    b.label("lblK", 150, 84, 15, "k:")
    b.dropdown("k", 165, 84, 85, K_VALUES)
    b.check("standardize", 22, 98, 228, "Standardise variables before computing distances",
            opts.standardize)
    b.label("lblMiceMethod", 22, 130, 60, "Method:")
    b.dropdown("miceMethod", 85, 130, 165, MICE_METHODS)
    b.label("lblM", 22, 146, 60, "Imputations:")
    b.dropdown("m", 85, 146, 45, M_VALUES)
    b.label("lblIter", 140, 146, 50, "Iterations:")
    b.dropdown("iterations", 195, 146, 55, ITERATIONS)
    b.label("lblSeed", 22, 162, 60, "Random seed:")
    b.edit("seed", 85, 162, 45, str(opts.seed))
    b.label("lblSeedHint", 135, 162, 115, "(empty = different each run)")

    b.line("lineOut", 180, "Output")
    b.radio("newSheet", 10, 192, 240, "Completed data on a new sheet", True)
    b.radio("inPlace", 10, 204, 240, "Fill the missing cells in the input range (not for MICE)",
            False)
    b.check("highlight", 10, 216, 240, "Highlight imputed cells", opts.highlight)
    b.ok_cancel(234)

    dialog = b.create(parent)
    d.select(dialog, "k", K_VALUES, opts.k)
    d.select(dialog, "miceMethod", MICE_METHODS, opts.mice_method)
    d.select(dialog, "m", M_VALUES, opts.m)
    d.select(dialog, "iterations", ITERATIONS, opts.iterations)

    def validate(result):
        if result.seed == "invalid":
            return "The random seed must be a whole number (or empty)."
        if result.method == "mice" and result.in_place:
            return ("Multiple imputation creates several datasets, so it always writes to a "
                    "new sheet. Choose 'Completed data on a new sheet'.")
        return None

    return d.run(ctx, doc, parent, dialog, IMPUTE_TITLE, _read_impute, validate)


def _read_impute(dialog):
    seed_text = dialog.getControl("seed").getText().strip()
    try:
        seed = int(seed_text) if seed_text else None
    except ValueError:
        seed = "invalid"
    method = next(m for m in METHODS if d.checked(dialog, m))
    return ImputeOptions(
        range_text=dialog.getControl("range").getText(),
        has_labels=d.checked(dialog, "labels"),
        method=method,
        k=d.picked(dialog, "k", K_VALUES),
        standardize=d.checked(dialog, "standardize"),
        mice_method=d.picked(dialog, "miceMethod", MICE_METHODS),
        m=d.picked(dialog, "m", M_VALUES),
        iterations=d.picked(dialog, "iterations", ITERATIONS),
        seed=seed,
        in_place=d.checked(dialog, "inPlace"),
        highlight=d.checked(dialog, "highlight"))
