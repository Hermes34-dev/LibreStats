"""Dialogs for Data Visualization."""

from . import calc
from . import dialogs as d
from .prep_dialog import SCORES
from .viz_calc import DEFAULT_COLOR, HistogramOptions, QQOptions

COLORS = ((DEFAULT_COLOR, "Blue"), (0x729FCF, "Light blue"), (0x579D1C, "Green"),
          (0xFF8000, "Orange"), (0xC9211E, "Red"), (0x7E0021, "Burgundy"),
          (0x808080, "Grey"), (0x000000, "Black"))


def show_histogram(ctx, doc, parent):
    """Show the dialog; returns HistogramOptions, or None if cancelled."""
    title = "LibreStats - Histogram"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 300)
    b.range_and_labels(range_text, labels)

    b.line("lineBins", 42, "Bins (all selected cells are pooled; text is ignored)")
    b.radio("binsAuto", 10, 54, 240, "Automatic (Sturges' rule: ceil(log2 n) + 1 bins)", True)
    b.radio("binsCount", 10, 68, 120, "Number of bins:", False)
    b.radio("binsWidth", 10, 82, 120, "Bin width:", False)
    b.edit("k", 135, 68, 40, "10")
    b.edit("width", 135, 82, 60, "")

    b.line("lineRange", 100, "Range (leave empty to use the data's minimum and maximum)")
    b.label("lblMin", 10, 112, 45, "Minimum:")
    b.edit("min", 55, 112, 60, "")
    b.label("lblMax", 135, 112, 45, "Maximum:")
    b.edit("max", 180, 112, 60, "")

    b.line("lineY", 130, "Y axis")
    b.radio("count", 10, 142, 120, "Count", True)
    b.radio("percent", 135, 142, 115, "Percent", False)

    b.line("lineChart", 160, "Chart (leave titles empty for automatic ones)")
    b.label("lblTitle", 10, 172, 55, "Title:")
    b.edit("title", 70, 172, 180, "")
    b.label("lblXTitle", 10, 188, 55, "X axis title:")
    b.edit("xTitle", 70, 188, 180, "")
    b.label("lblYTitle", 10, 204, 55, "Y axis title:")
    b.edit("yTitle", 70, 204, 180, "")
    b.label("lblColor", 10, 220, 55, "Bar colour:")
    b.dropdown("color", 70, 220, 80, COLORS)

    b.line("lineOut", 238, "Output: frequency table and chart")
    b.radio("outNew", 10, 250, 240, "New sheet", True)
    b.radio("outBeside", 10, 264, 240, "To the right of the input range", False)
    b.ok_cancel(282)

    dialog = b.create(parent)
    d.select(dialog, "color", COLORS, DEFAULT_COLOR)
    problems = {}

    def read(dlg):
        text = lambda name: dlg.getControl(name).getText()  # noqa: E731
        opts = HistogramOptions(
            range_text=text("range"), has_labels=d.checked(dlg, "labels"),
            bins={"binsAuto": "auto", "binsCount": "count", "binsWidth": "width"}[
                d.radio_choice(dlg, ("binsAuto", "binsCount", "binsWidth"))],
            percent=d.checked(dlg, "percent"),
            title=text("title"), x_title=text("xTitle"), y_title=text("yTitle"),
            color=d.picked(dlg, "color", COLORS), beside=d.checked(dlg, "outBeside"))
        problems.pop("msg", None)
        try:
            if opts.bins == "count":
                k = d.number(text("k"), "The number of bins")
                if k is None or k != int(k) or k < 1:
                    raise ValueError("The number of bins must be a whole number of 1 or more.")
                opts.k = int(k)
            if opts.bins == "width":
                opts.width = d.number(text("width"), "The bin width")
                if opts.width is None or opts.width <= 0:
                    raise ValueError("Enter a positive bin width.")
            opts.min = d.number(text("min"), "The minimum")
            opts.max = d.number(text("max"), "The maximum")
            if opts.min is not None and opts.max is not None and opts.min >= opts.max:
                raise ValueError("The minimum must be smaller than the maximum.")
        except ValueError as exc:
            problems["msg"] = str(exc)
        return opts

    return d.run(ctx, doc, parent, dialog, title, read, lambda _: problems.get("msg"))


def show_qq(ctx, doc, parent):
    """Show the dialog; returns QQOptions, or None if cancelled."""
    title = "LibreStats - Normal Q-Q Plot"
    range_text, labels = d.initial_range(doc)
    b = d.Builder(ctx, title, 260, 260)
    b.range_and_labels(range_text, labels)

    b.line("lineGroups", 42, "Groups")
    b.label("lblGroups", 10, 54, 85, "Group labels (optional):")
    b.edit("groupRange", 100, 54, 150, "")
    b.label("lblGroupsHint", 10, 68, 240, "Empty: one plot per column. Else: one column of "
            "values, per group.")

    b.line("linePlot", 86, "Plot")
    b.label("lblPositions", 10, 98, 85, "Plotting positions:")
    b.dropdown("positions", 100, 98, 150, SCORES)
    b.label("lblLine", 10, 116, 85, "Reference line:")
    b.radio("quartiles", 100, 116, 150, "through the quartiles (as R qqline)", True)
    b.radio("meansd", 100, 130, 150, "mean and SD", False)
    b.radio("none", 100, 144, 150, "none", False)

    b.line("lineChart", 162, "Chart (leave titles empty for automatic ones)")
    b.label("lblTitle", 10, 174, 55, "Title:")
    b.edit("title", 70, 174, 180, "")
    b.label("lblXTitle", 10, 190, 55, "X axis title:")
    b.edit("xTitle", 70, 190, 180, "")
    b.label("lblYTitle", 10, 206, 55, "Y axis title:")
    b.edit("yTitle", 70, 206, 180, "")
    b.label("lblColor", 10, 222, 55, "Point colour:")
    b.dropdown("color", 70, 222, 80, COLORS)
    b.ok_cancel(242)

    dialog = b.create(parent)
    d.select(dialog, "positions", SCORES, "blom")
    d.select(dialog, "color", COLORS, DEFAULT_COLOR)
    problems = {}

    def read(dlg):
        text = lambda name: dlg.getControl(name).getText()  # noqa: E731
        opts = QQOptions(
            range_text=text("range"), has_labels=d.checked(dlg, "labels"),
            group_range=text("groupRange"), positions=d.picked(dlg, "positions", SCORES),
            line=d.radio_choice(dlg, ("quartiles", "meansd", "none")),
            title=text("title"), x_title=text("xTitle"), y_title=text("yTitle"),
            color=d.picked(dlg, "color", COLORS))
        problems.pop("msg", None)
        if opts.group_range.strip():
            try:
                calc.resolve_range(doc, opts.group_range)
            except calc.RangeError as exc:
                problems["msg"] = "Group labels: " + str(exc)
        return opts

    return d.run(ctx, doc, parent, dialog, title, read, lambda _: problems.get("msg"))


DIALOGS = {"histogram": show_histogram, "qq": show_qq}
