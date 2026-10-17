"""Shared helpers for building LibreStats dialogs in code (sizes in dialog units)."""

from . import calc

PUSH_OK, PUSH_CANCEL = 1, 2  # com.sun.star.awt.PushButtonType


class Builder:
    """Adds controls to a dialog model. Radio buttons inserted one after another
    form a group; any other control inserted in between starts a new group."""

    def __init__(self, ctx, title, width, height):
        self.ctx = ctx
        self.model = ctx.ServiceManager.createInstanceWithContext(
            "com.sun.star.awt.UnoControlDialogModel", ctx)
        self.model.Title, self.model.Width, self.model.Height = title, width, height

    def _add(self, kind, name, x, y, w, h, **props):
        ctl = self.model.createInstance(f"com.sun.star.awt.UnoControl{kind}Model")
        ctl.PositionX, ctl.PositionY, ctl.Width, ctl.Height = x, y, w, h
        for key, value in props.items():
            setattr(ctl, key, value)
        self.model.insertByName(name, ctl)

    def label(self, name, x, y, w, text):
        self._add("FixedText", name, x, y + 2, w, 10, Label=text)

    def line(self, name, y, text):
        self._add("FixedLine", name, 6, y, self.model.Width - 12, 10, Label=text)

    def check(self, name, x, y, w, text, state):
        self._add("CheckBox", name, x, y, w, 12, Label=text, State=int(state))

    def radio(self, name, x, y, w, text, state):
        self._add("RadioButton", name, x, y, w, 12, Label=text, State=int(state))

    def edit(self, name, x, y, w, text):
        self._add("Edit", name, x, y, w, 12, Text=text)

    def dropdown(self, name, x, y, w, items):
        self._add("ListBox", name, x, y, w, 12, Dropdown=True,
                  StringItemList=tuple(label for _, label in items))

    def range_and_labels(self, range_text, has_labels):
        """Input range field and 'first row contains names' box at the top."""
        self.label("lblRange", 10, 8, 55, "Input range:")
        self.edit("range", 70, 8, self.model.Width - 80, range_text)
        self.check("labels", 10, 24, 240, "First row contains variable names", has_labels)

    def multiline(self, name, x, y, w, h, text):
        self._add("Edit", name, x, y, w, h, Text=text, MultiLine=True, VScroll=True)

    def button(self, name, x, y, w, text):
        self._add("Button", name, x, y, w, 14, Label=text)

    def output_choice(self, y):
        """'New sheet' / 'right of the input range' radio pair (12 + 30 units high)."""
        self.line("lineOut", y, "Output")
        self.radio("outNew", 10, y + 12, 240, "New sheet: copy of the data plus the new columns",
                   True)
        self.radio("outBeside", 10, y + 26, 240,
                   "New columns directly to the right of the input range", False)

    def ok_cancel(self, y):
        w = self.model.Width
        self._add("Button", "cancel", w - 116, y, 50, 14, Label="Cancel", PushButtonType=PUSH_CANCEL)
        self._add("Button", "ok", w - 60, y, 50, 14, Label="OK", PushButtonType=PUSH_OK,
                  DefaultButton=True)

    def create(self, parent):
        smgr = self.ctx.ServiceManager
        dialog = smgr.createInstanceWithContext("com.sun.star.awt.UnoControlDialog", self.ctx)
        dialog.setModel(self.model)
        toolkit = smgr.createInstanceWithContext("com.sun.star.awt.Toolkit", self.ctx)
        dialog.createPeer(toolkit, parent)
        return dialog


def index(items, value):
    return next(i for i, (v, _) in enumerate(items) if v == value)


def select(dialog, name, items, value):
    dialog.getControl(name).selectItemPos(index(items, value), True)


def checked(dialog, name):
    return dialog.getControl(name).getState() == 1


def picked(dialog, name, items):
    pos = dialog.getControl(name).getSelectedItemPos()
    return items[pos][0] if pos >= 0 else items[0][0]


def number(text, what):
    """Parse a number typed in a dialog (accepts a decimal comma); None if empty."""
    text = text.strip()
    if not text:
        return None
    try:
        return float(text.replace(",", "."))
    except ValueError:
        raise ValueError(f"{what} must be a number (got '{text}').") from None


def radio_choice(dialog, names):
    return next(n for n in names if checked(dialog, n))


def initial_range(doc):
    """Range text from the current selection and a guess whether row 1 holds names."""
    text = calc.selection_text(doc)
    try:
        first = calc.resolve_range(doc, text).getDataArray()[0]
    except calc.RangeError:
        return text, True
    return text, any(isinstance(v, str) and v.strip() for v in first)


def run(ctx, doc, parent, dialog, title, read, validate=None):
    """Execute until the input is valid; returns read(dialog), or None if cancelled."""
    try:
        while True:
            if dialog.execute() != 1:
                return None
            result = read(dialog)
            problem = validate(result) if validate else None
            if problem is None:
                try:
                    calc.resolve_range(doc, result.range_text)
                except calc.RangeError as exc:
                    problem = str(exc)
            if problem is None:
                return result
            calc.message_box(ctx, parent, title, problem)
    finally:
        dialog.dispose()
