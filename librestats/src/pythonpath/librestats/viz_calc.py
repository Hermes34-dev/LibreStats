"""Data visualisation: tables of plotted values and the Calc charts built from them."""

from dataclasses import dataclass

import uno
from com.sun.star.awt import Point, Rectangle

from . import assumptions_calc, calc
from . import histogram as h
from . import qq

DEFAULT_COLOR = 0x004586  # LibreOffice's default chart blue
LINE_COLOR = 0x333333      # reference lines
CHART_SIZE = (16000, 9000)  # 1/100 mm
QQ_SIZE = (12000, 10000)
LINE_SOLID = uno.Enum("com.sun.star.drawing.LineStyle", "SOLID")
LINE_NONE = uno.Enum("com.sun.star.drawing.LineStyle", "NONE")


@dataclass
class HistogramOptions:
    range_text: str = ""
    has_labels: bool = True
    bins: str = "auto"          # auto | count | width
    k: int = 10
    width: float = None
    min: float = None           # None: the data's own minimum
    max: float = None
    percent: bool = False
    title: str = ""             # empty titles are filled in automatically
    x_title: str = ""
    y_title: str = ""
    color: int = DEFAULT_COLOR
    beside: bool = False


@dataclass
class QQOptions:
    range_text: str = ""
    has_labels: bool = True
    group_range: str = ""
    positions: str = "blom"     # blom | vdw | tukey | rankit
    line: str = "quartiles"     # quartiles | meansd | none
    title: str = ""
    x_title: str = ""
    y_title: str = ""
    color: int = DEFAULT_COLOR


def _darker(color, factor=0.6):
    r, g, b = (color >> 16) & 255, (color >> 8) & 255, color & 255
    return (int(r * factor) << 16) | (int(g * factor) << 8) | int(b * factor)


def _write_block(sheet, col, row, rows):
    """Write rows (padded to equal width) with the top-left cell at (col, row)."""
    width = max(len(r) for r in rows)
    data = tuple(tuple(calc._cell_value(v) for v in list(r) + [""] * (width - len(r)))
                 for r in rows)
    block = sheet.getCellRangeByPosition(col, row, col + width - 1, row + len(rows) - 1)
    block.setDataArray(data)
    return block


def _add_chart(doc, sheet, name, table, pos, size, row_headers, title, x_title, y_title):
    """Insert a chart of table at pos (1/100 mm) without legend; return the chart model."""
    # chart names must be unique in the whole document, not only on the sheet
    sheets = doc.getSheets()
    taken = {n for s in range(sheets.getCount())
             for n in sheets.getByIndex(s).getCharts().getElementNames()}
    base, i = name, 1
    while name in taken:
        i += 1
        name = f"{base} {i}"
    charts = sheet.getCharts()
    charts.addNewByName(name, Rectangle(pos.X, pos.Y, *size), (table.getRangeAddress(),),
                        True, row_headers)
    chart = charts.getByName(name).getEmbeddedObject()
    chart.HasMainTitle = True
    chart.getTitle().String = title
    chart.HasLegend = False
    return chart


def _axis_titles(chart, x_title, y_title):
    diagram = chart.getDiagram()
    diagram.HasXAxisTitle = True
    diagram.getXAxisTitle().String = x_title
    diagram.HasYAxisTitle = True
    diagram.getYAxisTitle().String = y_title


def _series(chart):
    chart_type = chart.getFirstDiagram().getCoordinateSystems()[0].getChartTypes()[0]
    return chart_type, chart_type.getDataSeries()


def insert_column_chart(doc, sheet, name, table, pos, title, x_title, y_title, color, gap=0):
    """Column chart of a two-column table (labels, values; header row)."""
    chart = _add_chart(doc, sheet, name, table, pos, CHART_SIZE, True, title, x_title, y_title)
    _axis_titles(chart, x_title, y_title)
    chart_type, series = _series(chart)
    chart_type.GapwidthSequence = (gap,)
    for s in series:
        s.Color = color
        s.BorderStyle = LINE_SOLID
        s.BorderColor = _darker(color)
    return chart


def insert_scatter_chart(doc, sheet, name, table, pos, size, title, x_title, y_title, color,
                         line_color=None):
    """XY chart of a table (x, y[, line]; header row): y as points, line as a line."""
    chart = _add_chart(doc, sheet, name, table, pos, size, False, title, x_title, y_title)
    chart.setDiagram(chart.createInstance("com.sun.star.chart.XYDiagram"))
    _axis_titles(chart, x_title, y_title)
    _, series = _series(chart)
    points = series[0]
    points.LineStyle = LINE_NONE
    symbol = points.Symbol
    symbol.Style = uno.Enum("com.sun.star.chart2.SymbolStyle", "STANDARD")
    symbol.StandardSymbol = 1   # circle
    symbol.FillColor = symbol.BorderColor = color
    points.Symbol = symbol
    points.Color = color
    for s in series[1:]:
        symbol = s.Symbol
        symbol.Style = uno.Enum("com.sun.star.chart2.SymbolStyle", "NONE")
        s.Symbol = symbol
        s.LineStyle = LINE_SOLID
        s.Color = LINE_COLOR if line_color is None else line_color
        s.LineWidth = 35
    return chart


def _variable_name(header):
    names = [str(v).strip() for v in header if str(v).strip()]
    return names[0] if len(names) == 1 else ""


def run_histogram(doc, opts):
    cell_range = calc.resolve_range(doc, opts.range_text)
    addr = cell_range.getRangeAddress()
    data = calc.used_data(cell_range)
    name = _variable_name(data[0]) if opts.has_labels and data else ""
    values, text = h.numeric_values(data[1:] if opts.has_labels else data)
    hist = h.histogram(values, opts.bins, opts.k, opts.width, opts.min, opts.max)

    y_label = "Percent" if opts.percent else "Count"
    freq = hist.percents() if opts.percent else hist.counts
    table = [["Bin", y_label]] + [[lab, f] for lab, f in zip(hist.labels(), freq)]

    info = [["Values counted", hist.n],
            ["Bins", hist.rule],
            ["Bin width", hist.width],
            ["Range", f"{h.fmt(hist.edges[0])} to {h.fmt(hist.edges[-1])}"],
            ["Intervals", "[from, to): the last bin includes its upper limit"]]
    if hist.excluded:
        info.append(["Outside the range (not counted)", hist.excluded])
    if text:
        info.append(["Text cells ignored", text])
    if opts.percent:
        info.append(["Percent of", f"the {hist.n} values counted"])
    info += [["Note", n] for n in hist.notes]

    title = opts.title.strip() or (f"Histogram of {name}" if name else "Histogram")
    x_title = opts.x_title.strip() or name or "Value"
    y_title = opts.y_title.strip() or y_label
    rows = table + [[]] + info

    with calc.undo_context(doc, "LibreStats: Histogram"):
        if opts.beside:
            sheet = cell_range.getSpreadsheet()
            col, row = addr.EndColumn + 2, addr.StartRow
            target = sheet.getCellRangeByPosition(col, row, col + 1, row + len(rows) - 1)
            if target.queryContentCells(1 | 2 | 4 | 16).getCount():
                raise calc.RangeError("The cells to the right of the input range are not "
                                      "empty. Clear them or choose 'New sheet'.")
        else:
            sheet = calc.new_sheet(doc, "Histogram")
            col, row = 0, 0
        table_range = _write_block(sheet, col, row, table)
        _write_block(sheet, col, row + len(table) + 1, info)
        sheet.getCellRangeByPosition(col, row, col + 1, row).CharWeight = calc.BOLD
        sheet.getCellRangeByPosition(col, row, col + 1, row).getColumns().OptimalWidth = True
        insert_column_chart(doc, sheet, "Histogram", table_range,
                            sheet.getCellByPosition(col + 3, row).Position, title, x_title,
                            y_title, opts.color)
    doc.getCurrentController().setActiveSheet(sheet)
    return sheet.getName()


LINE_NAMES = {"quartiles": "through the 1st and 3rd quartiles",
              "meansd": "intercept = mean, slope = SD", "none": "none"}


def run_qq(doc, opts):
    samples = assumptions_calc.read_samples(doc, opts.range_text, opts.has_labels,
                                            opts.group_range)
    plots, skipped = [], []
    for s in samples.items:
        try:
            plots.append((s, qq.qq_normal(s.values, opts.positions, opts.line)))
        except qq.QQError as exc:
            skipped.append(f"{s.name}: {exc}")
    if not plots:
        raise calc.RangeError("No Q-Q plot could be drawn. " + "; ".join(skipped))

    with calc.undo_context(doc, "LibreStats: Q-Q plot"):
        sheet = calc.new_sheet(doc, "Q-Q plot")
        width = 2 if opts.line == "none" else 3
        tables, col = [], 0
        for s, res in plots:
            info = [[s.name], ["n", len(res.sample)],
                    ["Plotting positions", qq.describe(opts.positions)],
                    ["Reference line", LINE_NAMES[opts.line]]]
            if res.slope is not None:
                info += [["Intercept", res.intercept], ["Slope", res.slope]]
            info += [["Note", n] for n in res.notes]
            _write_block(sheet, col, 0, info)
            sheet.getCellByPosition(col, 0).CharWeight = calc.BOLD
            top = len(info) + 1
            header = ["Normal quantile", s.name] + (["Reference line"] if width == 3 else [])
            line = res.line()
            rows = [header] + [[z, y] + ([line[i]] if line else [])
                               for i, (z, y) in enumerate(zip(res.theoretical, res.sample))]
            tables.append(_write_block(sheet, col, top, rows))
            sheet.getCellRangeByPosition(col, top, col + width - 1, top).CharWeight = calc.BOLD
            sheet.getCellRangeByPosition(col, 0, col + width - 1, 0).getColumns() \
                .OptimalWidth = True
            col += width + 1

        origin = sheet.getCellByPosition(col, 0).Position
        for i, ((s, _), table) in enumerate(zip(plots, tables)):
            if opts.title.strip():
                title = opts.title.strip() + (f" ({s.name})" if len(plots) > 1 else "")
            else:
                title = f"Normal Q-Q plot of {s.name}"
            y_title = opts.y_title.strip() or samples.variable or s.name
            pos = Point(origin.X, origin.Y + i * (QQ_SIZE[1] + 500))
            insert_scatter_chart(doc, sheet, "Q-Q plot", table, pos, QQ_SIZE, title,
                                 opts.x_title.strip() or "Theoretical quantiles (standard normal)",
                                 y_title, opts.color)
        if skipped:
            note_row = max(t.getRangeAddress().EndRow for t in tables) + 2
            _write_block(sheet, 0, note_row, [["Not plotted: " + m] for m in skipped])
    doc.getCurrentController().setActiveSheet(sheet)
    return sheet.getName()


RUNNERS = {"histogram": run_histogram, "qq": run_qq}
