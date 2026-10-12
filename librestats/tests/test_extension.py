"""Headless integration test of the installed LibreStats extension.

Run with LibreOffice's Python after installing the extension:
    /opt/libreoffice26.8/program/python tests/test_extension.py
"""

import csv
import os
import random
import signal
import subprocess
import sys
import time

import uno
from com.sun.star.beans import NamedValue, PropertyValue

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src",
                                "pythonpath"))
from librestats import compositional, transform  # noqa: E402  (pure Python, for expected values)

SOFFICE = "/opt/libreoffice26.8/program/soffice"
PIPE = "librestats_test"
HIGHLIGHT = 0xFFCC99
IMPUTED = 0xCCE5FF
AIRQUALITY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "airquality.csv")

failures = []


def check(condition, message):
    print(("PASS  " if condition else "FAIL  ") + message)
    if not condition:
        failures.append(message)


def connect(timeout=60):
    local = uno.getComponentContext()
    resolver = local.ServiceManager.createInstanceWithContext(
        "com.sun.star.bridge.UnoUrlResolver", local)
    deadline = time.time() + timeout
    while True:
        try:
            return resolver.resolve(f"uno:pipe,name={PIPE};urp;StarOffice.ComponentContext")
        except Exception:
            if time.time() > deadline:
                raise
            time.sleep(0.5)


def fill_data(sheet):
    """A: Score (outliers in rows 9 and 13, plus a text and an empty cell);
    B, C: X and Y on a line, row 22 far off the line but within both ranges."""
    rnd = random.Random(7)
    rows = [("Score", "X", "Y")]
    for i in range(1, 41):
        score = round(rnd.gauss(50, 5), 2)
        x = float(i)
        y = round(2 * x + rnd.gauss(0, 1), 2)
        rows.append((score, x, y))
    rows[8] = (120.0,) + rows[8][1:]    # A9: high outlier
    rows[12] = (-20.0,) + rows[12][1:]  # A13: low outlier
    rows[21] = (rows[21][0], 21.0, 5.0)  # row 22: Y should be ~42
    rows[30] = ("n/a",) + rows[30][1:]  # A31: text, ignored
    rows[31] = ("",) + rows[31][1:]     # A32: empty, ignored
    sheet.getCellRangeByName("A1:C41").setDataArray(tuple(tuple(r) for r in rows))


def report_text(sheet):
    cursor = sheet.createCursor()
    cursor.gotoEndOfUsedArea(False)
    used = sheet.getCellRangeByPosition(0, 0, cursor.RangeAddress.EndColumn,
                                        cursor.RangeAddress.EndRow)
    return used.getDataArray()


def load_airquality(sheet):
    """airquality (R) in A1:F154; missing = empty, except A6 holds the text 'NA'."""
    with open(AIRQUALITY) as f:
        rows = [r[1:] for r in csv.reader(f)]
    data = [tuple(rows[0])] + [tuple(float(v) if v else "" for v in r) for r in rows[1:]]
    data[5] = ("NA",) + data[5][1:]
    sheet.getCellRangeByName("A1:F154").setDataArray(tuple(data))
    return rows


def find_row(rows, first):
    return next((r for r in rows if r and r[0] == first), None)


def check_missing(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    load_airquality(data)
    base = (NamedValue("Document", doc), NamedValue("Range", "$Sheet1.$A$1:$F$154"))

    name = handler.execute(base + (NamedValue("Command", "missing_diagnose"),))
    rows = report_text(doc.Sheets.getByName(name))
    check(name == "Missing data", f"Diagnostics sheet created ({name!r})")
    check(find_row(rows, "Text cells counted as missing")[1] == 1.0, "Text 'NA' counted as missing")
    check(find_row(rows, "Ozone")[2] == 37.0, "Ozone: 37 missing")
    header = rows.index(next(r for r in rows if r[0] == "Chi-square"))
    stat, df, p, patterns = rows[header + 1][:4]
    check(abs(stat - 35.10613) < 1e-3 and df == 14 and abs(p - 0.0014178) < 1e-6
          and patterns == 4, f"Little's test matches naniar (chi2={stat:.5f}, df={df:g}, p={p:.6f})")

    name = handler.execute(base + (NamedValue("Command", "impute"), NamedValue("Method", "mean")))
    imputed = doc.Sheets.getByName(name)
    ozone = [float(r[0]) for r in load_rows() if r[0]]
    check(abs(imputed.getCellRangeByName("A6").getValue() - sum(ozone) / len(ozone)) < 1e-9,
          "Mean imputation: A6 = mean Ozone")
    check(imputed.getCellRangeByName("A6").CellBackColor == IMPUTED, "Imputed cell highlighted")
    check(imputed.getCellRangeByName("A2").CellBackColor != IMPUTED, "Observed cell not highlighted")
    check(imputed.getCellRangeByName("A1").getString() == "Ozone", "Header row copied")

    name = handler.execute(base + (NamedValue("Command", "impute"), NamedValue("Method", "mice"),
                                   NamedValue("M", 3), NamedValue("Iterations", 5),
                                   NamedValue("Seed", 1)))
    mice_rows = report_text(doc.Sheets.getByName(name))
    block = [r[:8] for r in mice_rows[1:1 + 3 * 153]]
    check(len([r for r in block if r[0] in (1.0, 2.0, 3.0)]) == 3 * 153, "MICE: 3 stacked datasets")
    check(all(isinstance(v, float) for r in block for v in r), "MICE: no missing values left")
    observed = {float(v) for v in ozone}
    check(block[4][2] in observed and block[4][1] == 6.0, "MICE (PMM): row 6 Ozone is an observed value")
    check(doc.Sheets.getByName(name).getCellRangeByName("C6").CellBackColor == IMPUTED,
          "MICE: imputed cell highlighted")

    name = handler.execute(base + (NamedValue("Command", "impute"), NamedValue("Method", "knn"),
                                   NamedValue("InPlace", True)))
    values = data.getCellRangeByName("A2:F154").getDataArray()
    check(all(isinstance(v, float) for r in values for v in r), "kNN in place: all cells filled")
    check(data.getCellRangeByName("A6").CellBackColor == IMPUTED, "kNN in place: A6 highlighted")
    check(name == "Imputation", "kNN in place: summary sheet created")

    sheets_before = doc.Sheets.getCount()
    for command in ("missingDiagnose", "missingImpute"):
        handler.trigger(command)
    check(doc.Sheets.getCount() == sheets_before, "Missing-data menu commands (cancelled dialogs) ran")
    doc.close(True)


PREP_DATA = [
    ("Conc", "Delta", "Group", "FA1", "FA2", "FA3"),
    (0.0, -2.5, "a", 20.0, 30.0, 50.0),
    (1.2, 0.4, "b", 25.0, 25.0, 50.0),
    (3.5, 1.1, "", 10.0, 60.0, 30.0),
    (2.2, -0.7, "c", 33.0, 33.0, 34.0),
    (9.8, 3.9, "a", 5.0, 15.0, 80.0),
    (4.1, 0.0, "b", 40.0, 40.0, 20.0),
    (0.6, -1.4, "a", 12.0, 18.0, 70.0),
    (5.5, 2.2, "c", 22.0, 28.0, 50.0),
]


def column_values(sheet, col, rows):
    return [sheet.getCellByPosition(col, r).getValue() for r in rows]


def prep(handler, doc, command, rng, **settings):
    args = (NamedValue("Command", command), NamedValue("Document", doc), NamedValue("Range", rng))
    return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))


def check_prep(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    data.getCellRangeByName("A1:F9").setDataArray(tuple(PREP_DATA))
    conc = [r[0] for r in PREP_DATA[1:]]
    delta = [r[1] for r in PREP_DATA[1:]]
    rows = range(1, 9)

    name = prep(handler, doc, "transform", "$Sheet1.$A$1:$A$9", method="log", offset="one")
    sheet = doc.Sheets.getByName(name)
    check(sheet.getCellRangeByName("B1").getString() == "ln(Conc+1)", "Log: header ln(Conc+1)")
    check(all(abs(a - b) < 1e-12 for a, b in zip(column_values(sheet, 1, rows),
                                                 [__import__("math").log(x + 1) for x in conc])),
          "Log: values ln(x+1)")

    name = prep(handler, doc, "transform", "$Sheet1.$B$1:$B$9", method="yeojohnson")
    lam = dict(transform.yeojohnson(delta).params)["Lambda (ML)"]
    rows_text = report_text(doc.Sheets.getByName(name))
    check(any(isinstance(v, float) and abs(v - lam) < 1e-9 for r in rows_text for v in r),
          f"Yeo-Johnson: lambda {lam:.4f} in report")

    try:
        prep(handler, doc, "transform", "$Sheet1.$A$1:$A$9", method="boxcox")
        check(False, "Box-Cox on data with a zero is refused")
    except Exception as exc:
        check("positive" in str(exc), "Box-Cox on data with a zero is refused")

    try:
        prep(handler, doc, "scale", "$Sheet1.$B$1:$B$9", method="z", Beside=True)
        check(False, "Writing beside the input refuses to overwrite column C")
    except Exception as exc:
        check("not empty" in str(exc) and data.getCellRangeByName("C1").getString() == "Group",
              "Writing beside the input refuses to overwrite column C")

    data.getCellRangeByName("H1:H9").setDataArray(tuple((v,) for v in ["Delta"] + delta))
    prep(handler, doc, "scale", "$Sheet1.$H$1:$H$9", method="pareto", Beside=True)
    expected = transform.scale(delta, "pareto").values
    check(data.getCellRangeByName("I1").getString() == "pareto(Delta)" and
          all(abs(a - b) < 1e-12 for a, b in zip(column_values(data, 8, rows), expected)),
          "Pareto scaling written beside the input")

    name = prep(handler, doc, "logratio", "$Sheet1.$D$1:$F$9", method="ilr")
    sheet = doc.Sheets.getByName(name)
    first = [r[3:6] for r in PREP_DATA[1:]][0]
    z1, z2 = compositional.ilr_pivot(list(first))
    check(abs(sheet.getCellRangeByName("D2").getValue() - z1) < 1e-12 and
          abs(sheet.getCellRangeByName("E2").getValue() - z2) < 1e-12 and
          sheet.getCellRangeByName("F2").getString() == "", "ILR: 2 pivot coordinates for 3 parts")

    name = prep(handler, doc, "dummies", "$Sheet1.$C$1:$C$9")
    sheet = doc.Sheets.getByName(name)
    headers = [sheet.getCellByPosition(c, 0).getString() for c in range(1, 4)]
    check(headers == ["Group=a", "Group=b", "Group=c"], f"Dummies: headers {headers}")
    check(sheet.getCellRangeByName("B2").getValue() == 1.0 and
          sheet.getCellRangeByName("B4").getString() == "", "Dummies: values and missing")

    name = prep(handler, doc, "recode", "$Sheet1.$C$1:$C$9", rules="a; b -> ab\nc -> other")
    sheet = doc.Sheets.getByName(name)
    check([sheet.getCellByPosition(1, r).getString() for r in (1, 2, 3, 4)] ==
          ["ab", "ab", "", "other"], "Recode: rules applied")

    name = prep(handler, doc, "bins", "$Sheet1.$A$1:$A$9", method="quantile", k=4)
    codes = column_values(doc.Sheets.getByName(name), 1, rows)
    check(sorted(codes) == [1, 1, 2, 2, 3, 3, 4, 4], f"Quantile bins: 2 per bin ({codes})")

    name = prep(handler, doc, "ranks", "$Sheet1.$A$1:$A$9", scores=True)
    sheet = doc.Sheets.getByName(name)
    check(column_values(sheet, 1, rows) == transform.ranks(conc).values and
          all(abs(a - b) < 1e-12 for a, b in zip(column_values(sheet, 2, rows),
                                                 transform.normal_scores(conc).values)),
          "Ranks and Blom normal scores")

    name = prep(handler, doc, "presence", "$Sheet1.$A$1:$A$9", threshold=1.0)
    check(column_values(doc.Sheets.getByName(name), 1, rows) == [0, 1, 1, 1, 1, 1, 0, 1],
          "Presence/absence above 1")

    sheets_before = doc.Sheets.getCount()
    doc.CurrentController.setActiveSheet(data)
    doc.CurrentController.select(data.getCellRangeByName("C1:C9"))
    for command in ("transform", "scale", "logratio", "dummies", "presence", "bins", "recode",
                    "ranks"):
        handler.trigger(command)
    check(doc.Sheets.getCount() == sheets_before, "All 8 new dialogs open and cancel cleanly")
    doc.close(True)


def load_rows():
    with open(AIRQUALITY) as f:
        return [r[1:] for r in csv.reader(f)][1:]


def main():
    office = subprocess.Popen([SOFFICE, "--headless", "--invisible", "--nologo", "--norestore",
                               f"--accept=pipe,name={PIPE};urp;"], start_new_session=True)
    desktop = None
    try:
        ctx = connect()
        smgr = ctx.ServiceManager

        provider = smgr.createInstanceWithContext(
            "com.sun.star.configuration.ConfigurationProvider", ctx)
        menus = provider.createInstanceWithArguments(
            "com.sun.star.configuration.ConfigurationAccess",
            (PropertyValue(Name="nodepath",
                           Value="/org.openoffice.Office.Addons/AddonUI/OfficeMenuBar"),))
        check(menus.hasByName("org.librestats.menu"), "LibreStats menu registered")

        desktop = smgr.createInstanceWithContext("com.sun.star.frame.Desktop", ctx)
        doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                           (PropertyValue(Name="Hidden", Value=True),))
        data = doc.Sheets.getByIndex(0)
        fill_data(data)

        handler = smgr.createInstanceWithContext("org.librestats.Handler", ctx)
        check(handler is not None, "Handler service created")

        args = (NamedValue("Command", "outliers"), NamedValue("Document", doc),
                NamedValue("Range", "$Sheet1.$A$1:$C$41"), NamedValue("HasLabels", True),
                NamedValue("Methods", "tukey,mad,grubbs,dixon,mahalanobis"),
                NamedValue("GrubbsIterate", True), NamedValue("Highlight", True))
        sheet_name = handler.execute(args)
        check(sheet_name == "Outliers", f"Report sheet created ({sheet_name!r})")

        rows = report_text(doc.Sheets.getByName(sheet_name))
        flat = [str(v) for row in rows for v in row]
        for title in ("Tukey fences", "Robust z-scores", "Grubbs test", "Dixon test",
                      "Mahalanobis distance"):
            check(any(s.startswith(title) for s in flat), f"Report has section: {title}")
        check("A9" in flat and "A13" in flat, "Univariate outliers A9 and A13 reported")
        check(any(s.startswith("Dixon's test needs 3 to 30") for s in flat),
              "Dixon reports n > 30 limit for 38-value column")
        mahal_rows = [row for row in rows if row and row[0] == 22.0]
        check(bool(mahal_rows), "Mahalanobis flags row 22")
        score_n = [row for row in rows if row and row[0] == "Score" and row[1] == 38.0]
        check(bool(score_n), "Text and empty cells skipped (Score n = 38)")

        check(data.getCellRangeByName("A9").CellBackColor == HIGHLIGHT, "A9 highlighted")
        check(data.getCellRangeByName("C22").CellBackColor == HIGHLIGHT, "C22 highlighted")
        check(data.getCellRangeByName("A2").CellBackColor != HIGHLIGHT, "A2 not highlighted")

        undo = doc.getUndoManager()
        check(undo.getCurrentUndoActionTitle() == "LibreStats: Outlier analysis",
              "Analysis is a single undo step")
        undo.undo()
        # Headless, undo only reverts formatting (sheet inserts and cell contents
        # need a visible window), so check the report sheet's removal in the UI.
        check(data.getCellRangeByName("A9").CellBackColor != HIGHLIGHT, "Undo removes highlight")

        # Dialog path: headless office auto-cancels the dialog; must not hang or crash.
        sheets_before = doc.Sheets.getCount()
        doc.CurrentController.select(data.getCellRangeByName("A1"))
        handler.trigger("outliers")
        check(doc.Sheets.getCount() == sheets_before, "Menu command (cancelled dialog) ran")

        doc.close(True)

        check_missing(desktop, handler)
        check_prep(desktop, handler)
    finally:
        try:
            desktop.terminate()
            office.wait(timeout=30)
        except Exception:
            os.killpg(office.pid, signal.SIGKILL)  # soffice wrapper and soffice.bin

    print(f"\n{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
