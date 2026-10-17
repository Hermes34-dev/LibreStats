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
from librestats import assumptions, compositional, qq, transform  # noqa: E402  (pure Python, for expected values)

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


def chart_parts(sheet):
    chart = sheet.Charts.getByIndex(0).EmbeddedObject
    chart_type = chart.getFirstDiagram().getCoordinateSystems()[0].getChartTypes()[0]
    return chart, chart_type


def check_histogram(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    # 0..10 with text and an empty cell; 2, 4, 6, 8 sit on bin edges, 10 is the maximum
    rows = [("Length",)] + [(float(i),) for i in range(11)] + [("n/a",), ("",)]
    data.getCellRangeByName("A1:A14").setDataArray(tuple(rows))

    def hist(rng, **settings):
        args = (NamedValue("Command", "histogram"), NamedValue("Document", doc),
                NamedValue("Range", rng))
        return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))

    name = hist("$Sheet1.$A$1:$A$14", Bins="count", K=5, Color=0xC9211E)
    sheet = doc.Sheets.getByName(name)
    table = sheet.getCellRangeByName("A1:B6").getDataArray()
    check(name == "Histogram", f"Histogram sheet created ({name!r})")
    check([r[0] for r in table] == ["Bin", "0–2", "2–4", "4–6", "6–8", "8–10"],
          f"Histogram: bin labels {[r[0] for r in table[1:]]}")
    check([r[1] for r in table[1:]] == [2.0, 2.0, 2.0, 2.0, 3.0],
          "Histogram: edge values go up a bin, maximum in last bin")
    info = report_text(sheet)
    check(find_row(info, "Text cells ignored")[1] == 1.0, "Histogram: text cell reported")
    check(sheet.Charts.getCount() == 1, "Histogram: chart inserted")
    chart, chart_type = chart_parts(sheet)
    series = chart_type.getDataSeries()
    check(chart.getTitle().String == "Histogram of Length" and
          chart.Diagram.getXAxisTitle().String == "Length" and
          chart.Diagram.getYAxisTitle().String == "Count", "Histogram: automatic titles")
    check(tuple(chart_type.GapwidthSequence)[0] == 0, "Histogram: gap between bars is 0")
    check(len(series) == 1 and series[0].Color == 0xC9211E, "Histogram: one series, bar colour")
    check(doc.getUndoManager().getCurrentUndoActionTitle() == "LibreStats: Histogram",
          "Histogram is a single undo step")

    name = hist("$Sheet1.$A$1:$A$14", Bins="width", Width=3, Percent=True, Min=0, Max=10,
                Title="Lengths", XTitle="mm", YTitle="Share")
    sheet = doc.Sheets.getByName(name)
    table = sheet.getCellRangeByName("A1:B5").getDataArray()
    check(name == "Histogram 2" and table[0][1] == "Percent" and
          [r[0] for r in table[1:]] == ["0–3", "3–6", "6–9", "9–12"] and
          all(abs(r[1] - e) < 1e-12 for r, e in zip(table[1:], [300 / 11] * 3 + [200 / 11])),
          "Histogram: bin width 3 in percent")
    chart, _ = chart_parts(sheet)
    check((chart.getTitle().String, chart.Diagram.getXAxisTitle().String,
           chart.Diagram.getYAxisTitle().String) == ("Lengths", "mm", "Share"),
          "Histogram: custom titles")

    data.getCellRangeByName("D1:D2").setDataArray((("Single",), (4.5,)))
    hist("$Sheet1.$D$1:$D$2", Beside=True)
    check(data.getCellRangeByName("F2").getString() == "4.5–4.5" and
          data.getCellRangeByName("G2").getValue() == 1.0 and data.Charts.getCount() == 1,
          "Histogram: single value, one bin, written beside the data")
    try:
        hist("$Sheet1.$D$1:$D$2", Beside=True)
        check(False, "Histogram beside the data refuses to overwrite")
    except Exception as exc:
        check("not empty" in str(exc), "Histogram beside the data refuses to overwrite")

    for rng, what in (("$Sheet1.$J$1:$J$5", "empty selection"),
                      ("$Sheet1.$A$1", "text only")):
        try:
            hist(rng)
            check(False, f"Histogram of {what} is refused")
        except Exception as exc:
            check("no numeric values" in str(exc), f"Histogram of {what} is refused")

    sheets_before = doc.Sheets.getCount()
    handler.trigger("histogram")
    check(doc.Sheets.getCount() == sheets_before, "Histogram dialog opens and cancels cleanly")
    doc.close(True)


def check_assumptions(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    rnd = random.Random(11)
    normal = [round(rnd.gauss(50, 5), 2) for _ in range(30)]
    skewed = [round(rnd.expovariate(0.5), 3) for _ in range(30)]
    rows = [("Normal", "Skewed", "Group", "Value")]
    for i in range(30):  # D: Normal then Skewed values stacked, C: their group labels
        rows.append((normal[i], skewed[i], "", ""))
    rows[5] = ("n/a",) + rows[5][1:]  # A6: text, ignored
    stacked = [("ctrl", v) for v in normal[:15]] + [("treat", 3 * v) for v in normal[15:]]
    for i, (g, v) in enumerate(stacked):
        rows[i + 1] = rows[i + 1][:2] + (g, v)
    data.getCellRangeByName("A1:D31").setDataArray(tuple(rows))
    data.getCellRangeByName("F1:F3").setDataArray((("Tiny",), (1.0,), (2.0,)))
    normal_used = normal[:4] + normal[5:]

    def run(command, rng, **settings):
        args = (NamedValue("Command", command), NamedValue("Document", doc),
                NamedValue("Range", rng))
        return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))

    name = run("normality", "$Sheet1.$A$1:$B$31")
    rows = report_text(doc.Sheets.getByName(name))
    check(name == "Normality", f"Normality sheet created ({name!r})")
    for label, values in (("Normal", normal_used), ("Skewed", skewed)):
        ref = assumptions.shapiro_wilk(values)
        row = find_row(rows, label)
        check(row[1] == ref.n and abs(row[4] - ref.statistic) < 1e-12 and
              abs(row[5] - ref.p) < 1e-12, f"Shapiro-Wilk {label}: n = {ref.n}, W = "
              f"{ref.statistic:.4f}, p = {ref.p:.4g}")
    check(find_row(rows, "Skewed")[6].startswith("normality rejected"),
          "Shapiro-Wilk rejects normality of exponential data")
    check(find_row(rows, "Normal")[6] == "no evidence against normality",
          "Shapiro-Wilk keeps normality of normal data")
    check(any(str(r[0]).startswith("Note: 1 text cell") for r in rows if r),
          "Shapiro-Wilk: text cell reported")

    name = run("normality", "$Sheet1.$D$1:$D$31", GroupRange="$Sheet1.$C$1:$C$31")
    rows = report_text(doc.Sheets.getByName(name))
    ref = assumptions.shapiro_wilk([3 * v for v in normal[15:]])
    check(find_row(rows, "Group labels")[1] == "$Sheet1.$C$1:$C$31" and
          find_row(rows, "ctrl")[1] == 15 and abs(find_row(rows, "treat")[4] - ref.statistic)
          < 1e-12, "Shapiro-Wilk per group (values + group labels)")

    name = run("normality", "$Sheet1.$A$1:$F$31")
    rows = report_text(doc.Sheets.getByName(name))
    check(str(find_row(rows, "Tiny")[4]).startswith("not computed") and
          find_row(rows, "Normal")[1] == 29, "Shapiro-Wilk: too few values reported per column")

    name = run("variances", "$Sheet1.$A$1:$B$31")
    rows = report_text(doc.Sheets.getByName(name))
    ref = assumptions.levene([normal_used, skewed])
    header = rows.index(next(r for r in rows if r[0] == "F"))
    f, df1, df2, p = rows[header + 1][:4]
    check(name == "Equal variances" and abs(f - ref.statistic) < 1e-12 and (df1, df2) == (1, 57)
          and abs(p - ref.p) < 1e-12, f"Brown-Forsythe: F = {f:.4f}, df = {df1:g}, {df2:g}, "
          f"p = {p:.4g}")
    check("Brown-Forsythe" in rows[0][0], "Brown-Forsythe named in the report title")

    name = run("variances", "$Sheet1.$D$1:$D$31", GroupRange="$Sheet1.$C$1:$C$31",
               Center="mean")
    rows = report_text(doc.Sheets.getByName(name))
    ref = assumptions.levene([normal[:15], [3 * v for v in normal[15:]]], "mean")
    header = rows.index(next(r for r in rows if r[0] == "F"))
    check(abs(rows[header + 1][0] - ref.statistic) < 1e-12 and
          rows[header + 1][4].startswith("variances differ"),
          "Levene (mean) per group: SD x 3 gives unequal variances")

    for command, rng, settings, what in (
            ("variances", "$Sheet1.$A$1:$A$31", {}, "one group"),
            ("normality", "$Sheet1.$F$1:$F$3", {}, "too few values"),
            ("variances", "$Sheet1.$D$1:$D$31", {"GroupRange": "$Sheet1.$C$1:$C$20"},
             "group range of a different length")):
        try:
            run(command, rng, **settings)
            check(False, f"{command} refuses {what}")
        except Exception as exc:
            check(bool(str(exc)), f"{command} refuses {what}")

    sheets_before = doc.Sheets.getCount()
    for command in ("normality", "variances"):
        handler.trigger(command)
    check(doc.Sheets.getCount() == sheets_before, "Assumption dialogs open and cancel cleanly")
    doc.close(True)


def check_qq(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    rnd = random.Random(3)
    normal = [round(rnd.gauss(50, 5), 2) for _ in range(20)]
    skewed = [round(rnd.expovariate(0.5), 3) for _ in range(20)]
    rows = [("Normal", "Skewed", "Group", "Tiny")] + [(a, b, "g1" if i < 10 else "g2", "")
                                                    for i, (a, b) in enumerate(zip(normal, skewed))]
    rows[3] = ("n/a",) + rows[3][1:]          # A4: text, ignored
    rows[1] = rows[1][:3] + (1.0,)
    rows[2] = rows[2][:3] + (2.0,)
    data.getCellRangeByName("A1:D21").setDataArray(tuple(rows))
    normal_used = normal[:2] + normal[3:]

    def run(rng, **settings):
        args = (NamedValue("Command", "qq"), NamedValue("Document", doc), NamedValue("Range", rng))
        return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))

    name = run("$Sheet1.$A$1:$B$21", Color=0x579D1C)
    sheet = doc.Sheets.getByName(name)
    cells = report_text(sheet)
    header = next(i for i, r in enumerate(cells) if r[0] == "Normal quantile")
    ref = qq.qq_normal(normal_used)
    table = [r[:3] for r in cells[header + 1:header + 1 + len(ref.sample)]]
    check(name == "Q-Q plot" and cells[header][:3] == ("Normal quantile", "Normal", "Reference line"),
          f"Q-Q sheet created ({name!r}) with table header")
    check(all(abs(r[0] - z) < 1e-12 and r[1] == y and abs(r[2] - l) < 1e-12
              for r, z, y, l in zip(table, ref.theoretical, ref.sample, ref.line())) and
          len(table) == 19, "Q-Q: Blom quantiles, sorted values, quartile line (text skipped)")
    check(abs(find_row(cells, "Slope")[1] - ref.slope) < 1e-12, "Q-Q: slope reported")
    check(sheet.Charts.getCount() == 2, "Q-Q: one chart per column")
    chart, chart_type = chart_parts(sheet)
    series = chart_type.getDataSeries()
    check(chart_type.ChartType == "com.sun.star.chart2.ScatterChartType" and len(series) == 2,
          "Q-Q: XY scatter chart with points and line")
    check(series[0].LineStyle.value == "NONE" and series[0].Color == 0x579D1C and
          series[1].Symbol.Style.value == "NONE" and series[1].LineStyle.value == "SOLID",
          "Q-Q: points without lines, line without symbols, point colour")
    check(chart.getTitle().String == "Normal Q-Q plot of Normal" and
          chart.Diagram.getXAxisTitle().String == "Theoretical quantiles (standard normal)" and
          chart.Diagram.getYAxisTitle().String == "Normal", "Q-Q: automatic titles")

    name = run("$Sheet1.$B$1:$B$21", GroupRange="$Sheet1.$C$1:$C$21", Line="none",
               Positions="vdw", Title="Check")
    sheet = doc.Sheets.getByName(name)
    cells = report_text(sheet)
    header = next(i for i, r in enumerate(cells) if r[0] == "Normal quantile")
    ref = qq.qq_normal(skewed[:10], "vdw", "none")
    chart, chart_type = chart_parts(sheet)
    check(sheet.Charts.getCount() == 2 and cells[header][1] == "g1" and cells[header][2] == "" and
          abs(cells[header + 1][0] - ref.theoretical[0]) < 1e-12 and
          len(chart_type.getDataSeries()) == 1, "Q-Q per group, van der Waerden, no line")
    check(chart.getTitle().String == "Check (g1)" and
          chart.Diagram.getYAxisTitle().String == "Skewed", "Q-Q: custom title per group")

    name = run("$Sheet1.$A$1:$D$21")
    cells = report_text(doc.Sheets.getByName(name))
    notes = [str(r[0]) for r in cells if str(r[0]).startswith("Not plotted")]
    check(doc.Sheets.getByName(name).Charts.getCount() == 2 and
          notes == ["Not plotted: Group: A Q-Q plot needs at least 3 values.",
                    "Not plotted: Tiny: A Q-Q plot needs at least 3 values."],
          f"Q-Q: text column and column with 2 values skipped and reported ({notes})")
    try:
        run("$Sheet1.$D$1:$D$21")
        check(False, "Q-Q of 2 values is refused")
    except Exception as exc:
        check("at least 3" in str(exc), "Q-Q of 2 values is refused")

    sheets_before = doc.Sheets.getCount()
    handler.trigger("qq")
    check(doc.Sheets.getCount() == sheets_before, "Q-Q dialog opens and cancels cleanly")
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
        check_histogram(desktop, handler)
        check_qq(desktop, handler)
        check_assumptions(desktop, handler)
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
