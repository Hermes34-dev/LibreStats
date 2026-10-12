"""Helpers for working with Calc documents through UNO."""

import math
from contextlib import contextmanager

from com.sun.star.awt.MessageBoxType import ERRORBOX, INFOBOX
from com.sun.star.table import CellRangeAddress

BUTTONS_OK = 1  # com.sun.star.awt.MessageBoxButtons.BUTTONS_OK
BOLD = 150.0    # com.sun.star.awt.FontWeight.BOLD


class RangeError(ValueError):
    pass


def column_letter(index):
    """0 -> A, 25 -> Z, 26 -> AA."""
    letters = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        letters = chr(65 + rem) + letters
    return letters


def cell_name(column, row):
    return f"{column_letter(column)}{row + 1}"


def message_box(ctx, parent, title, text, error=True):
    toolkit = ctx.ServiceManager.createInstanceWithContext("com.sun.star.awt.Toolkit", ctx)
    box = toolkit.createMessageBox(parent, ERRORBOX if error else INFOBOX, BUTTONS_OK, title, text)
    box.execute()
    box.dispose()


def selection_text(doc):
    """Address of the current selection, expanded to the data region if it is one cell."""
    sel = doc.getCurrentSelection()
    if sel is None or not hasattr(sel, "getRangeAddress"):
        return ""
    addr = sel.getRangeAddress()
    if addr.StartColumn == addr.EndColumn and addr.StartRow == addr.EndRow:
        cursor = sel.getSpreadsheet().createCursorByRange(sel)
        cursor.collapseToCurrentRegion()
        return cursor.AbsoluteName
    return sel.AbsoluteName


def resolve_range(doc, text):
    """Turn a range like 'Sheet1.A1:C20' or 'A1:C20' (active sheet) into a cell range."""
    text = text.strip()
    if not text:
        raise RangeError("Please enter an input range.")
    conv = doc.createInstance("com.sun.star.table.CellRangeAddressConversion")
    active = doc.getCurrentController().getActiveSheet().getRangeAddress().Sheet
    conv.ReferenceSheet = active
    try:
        conv.UserInterfaceRepresentation = text
    except Exception:
        try:
            conv.PersistentRepresentation = text
        except Exception:
            raise RangeError(f"'{text}' is not a valid cell range.") from None
    a = conv.Address
    sheet = doc.getSheets().getByIndex(a.Sheet)
    return sheet.getCellRangeByPosition(a.StartColumn, a.StartRow, a.EndColumn, a.EndRow)


def new_sheet(doc, base_name):
    """Append a sheet named base_name (or 'base_name 2', ...) and return it."""
    sheets = doc.getSheets()
    name, i = base_name, 1
    while sheets.hasByName(name):
        i += 1
        name = f"{base_name} {i}"
    sheets.insertNewByName(name, sheets.getCount())
    return sheets.getByName(name)


def _cell_value(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, (int, float)):
        return float(v) if math.isfinite(v) else str(v)
    return str(v)


def write_rows(sheet, rows, bold_rows=()):
    """Write a list of rows (numbers / strings / None) starting at A1."""
    if not rows:
        return
    width = max(len(r) for r in rows) or 1
    data = tuple(tuple(_cell_value(v) for v in list(r) + [""] * (width - len(r)))
                 for r in rows)
    sheet.getCellRangeByPosition(0, 0, width - 1, len(rows) - 1).setDataArray(data)
    for r in bold_rows:
        sheet.getCellRangeByPosition(0, r, width - 1, r).CharWeight = BOLD
    sheet.getCellRangeByPosition(0, 0, width - 1, 0).getColumns().OptimalWidth = True


@contextmanager
def undo_context(doc, title):
    """Group all changes made inside the block into one undo step."""
    undo = doc.getUndoManager()
    undo.enterUndoContext(title)
    doc.lockControllers()
    try:
        yield
    finally:
        doc.unlockControllers()
        undo.leaveUndoContext()


def highlight_cells(doc, sheet, cells, color):
    """Set the background colour of (column, row) cells on sheet in one call."""
    if not cells:
        return
    ranges = doc.createInstance("com.sun.star.sheet.SheetCellRanges")
    index = sheet.getRangeAddress().Sheet
    for col, row in cells:
        ranges.addRangeAddress(CellRangeAddress(index, col, row, col, row), False)
    ranges.CellBackColor = color


def used_data(cell_range):
    """getDataArray() without trailing rows that are completely empty
    (e.g. when whole columns are selected)."""
    data = cell_range.getDataArray()
    end = len(data)
    while end and all(v == "" for v in data[end - 1]):
        end -= 1
    return data[:end]


class Report:
    """Rows for a report sheet; rows added with bold=True are formatted bold."""

    def __init__(self):
        self.rows, self.bold = [], []

    def add(self, *values, bold=False):
        if bold:
            self.bold.append(len(self.rows))
        self.rows.append(list(values))

    def blank(self):
        self.rows.append([])


def write_beside(sheet, data_rows, report):
    """Data table (header row bold) at A1 and the report one empty column to its right."""
    width = len(data_rows[0])
    merged = []
    for r in range(max(len(data_rows), len(report.rows))):
        left = data_rows[r] if r < len(data_rows) else [""] * width
        right = report.rows[r] if r < len(report.rows) else []
        merged.append(list(left) + [""] + list(right))
    write_rows(sheet, merged)
    sheet.getCellRangeByPosition(0, 0, width - 1, 0).CharWeight = BOLD
    for r in report.bold:
        sheet.getCellRangeByPosition(width + 1, r, width + 8, r).CharWeight = BOLD
