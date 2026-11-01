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
from librestats import assumptions, compare, compositional, posthoc, qq, transform  # noqa: E402  (pure Python, for expected values)

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


def check_compare(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    rnd = random.Random(23)
    pre = [round(rnd.gauss(50, 6), 2) for _ in range(20)]
    post = [round(v + rnd.gauss(2, 2), 2) for v in pre]
    late = [round(v + rnd.gauss(3, 4), 2) for v in pre]
    rows = [("Pre", "Post", "Late", "Group", "Value")]
    for i in range(20):  # E: values stacked by the labels in D (3 groups, unequal SDs)
        g = ("a", "b", "c")[i % 3]
        rows.append((pre[i], post[i], late[i], g, round(rnd.gauss(10 + i % 3, 1 + i % 3), 3)))
    rows[4] = (rows[4][0], "", rows[4][2]) + rows[4][3:]   # B5 empty: one incomplete row
    data.getCellRangeByName("A1:E21").setDataArray(tuple(rows))
    post_col = [r[1] if isinstance(r[1], float) else None for r in rows[1:]]
    groups = {}
    for r in rows[1:]:
        groups.setdefault(r[3], []).append(r[4])

    def run(command, rng, **settings):
        args = (NamedValue("Command", command), NamedValue("Document", doc),
                NamedValue("Range", rng))
        return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))

    def table_row(rows, header):
        return rows[rows.index(next(r for r in rows if r and r[0] == header)) + 1]

    name = run("ttestOne", "$Sheet1.$A$1:$B$21", Mu=50.0, Alternative="greater")
    rows = report_text(doc.Sheets.getByName(name))
    ref = compare.one_sample_t([v for v in post_col if v is not None], 50.0, "greater")
    row = find_row(rows, "Post")
    check(name == "One-sample t-test" and row[1] == 19 and abs(row[6] - ref.statistic) < 1e-12
          and abs(row[8] - ref.p) < 1e-12 and row[10] == "∞",
          f"One-sample t (greater): t = {ref.statistic:.4f}, p = {ref.p:.4g}, upper CI ∞")

    name = run("ttestWelch", "$Sheet1.$A$1:$B$21")
    rows = report_text(doc.Sheets.getByName(name))
    ref = compare.welch_t(pre, [v for v in post_col if v is not None])
    row = table_row(rows, "Difference")
    check(name == "Welch t-test" and abs(row[2] - ref.statistic) < 1e-12 and
          abs(row[3] - ref.df) < 1e-10 and abs(row[4] - ref.p) < 1e-12,
          f"Welch t: t = {ref.statistic:.4f}, df = {ref.df:.2f}, p = {ref.p:.4g}")

    name = run("ttestPaired", "$Sheet1.$A$1:$B$21")
    rows = report_text(doc.Sheets.getByName(name))
    ref = compare.paired_t(pre, post_col)
    row = table_row(rows, "Mean difference")
    check(name == "Paired t-test" and ref.n == 19 and abs(row[2] - ref.statistic) < 1e-12 and
          row[3] == 18 and abs(row[4] - ref.p) < 1e-12 and row[8].startswith("mean difference"),
          f"Paired t: t = {ref.statistic:.4f}, df = 18 (incomplete row dropped), p = {ref.p:.4g}")
    check(any(str(r[0]).startswith("Note: 1 row(s) with a missing value") for r in rows if r),
          "Paired t: incomplete row reported")

    name = run("anovaOneWay", "$Sheet1.$E$1:$E$21", GroupRange="$Sheet1.$D$1:$D$21")
    rows = report_text(doc.Sheets.getByName(name))
    ref = compare.one_way_anova(list(groups.values()))
    row = find_row(rows, "Between groups")
    check(name == "One-way ANOVA" and abs(row[4] - ref.statistic) < 1e-12 and row[2] == 2 and
          abs(row[5] - ref.p) < 1e-12 and find_row(rows, "Within groups")[2] == 17,
          f"One-way ANOVA per group: F = {ref.statistic:.4f}, df = 2, 17, p = {ref.p:.4g}")

    name = run("anovaWelch", "$Sheet1.$E$1:$E$21", GroupRange="$Sheet1.$D$1:$D$21")
    rows = report_text(doc.Sheets.getByName(name))
    ref = compare.welch_anova(list(groups.values()))
    row = table_row(rows, "F")
    check(name == "Welch ANOVA" and abs(row[0] - ref.statistic) < 1e-12 and
          abs(row[2] - ref.df2) < 1e-10 and abs(row[3] - ref.p) < 1e-12,
          f"Welch ANOVA: F = {ref.statistic:.4f}, df2 = {ref.df2:.2f}, p = {ref.p:.4g}")

    name = run("anovaRM", "$Sheet1.$A$1:$C$21")
    rows = report_text(doc.Sheets.getByName(name))
    ref = compare.rm_anova([[a, b, c] for a, b, c in zip(pre, post_col, late)])
    row = find_row(rows, "Conditions")
    gg = find_row(rows, "Greenhouse-Geisser")
    check(name == "RM ANOVA" and abs(row[4] - ref.statistic) < 1e-12 and
          find_row(rows, "Error")[2] == 36 and abs(gg[1] - ref.sphericity["gg"]) < 1e-12 and
          abs(gg[4] - ref.corrections[1][4]) < 1e-12,
          f"RM ANOVA: F = {ref.statistic:.4f}, GG ε = {ref.sphericity['gg']:.3f}, "
          f"p(GG) = {ref.corrections[1][4]:.4g}")
    check(abs(table_row(rows, "W")[0] - ref.sphericity["W"]) < 1e-12,
          "RM ANOVA: Mauchly's W reported")

    for command, rng, settings, what in (
            ("ttestWelch", "$Sheet1.$A$1:$C$21", {}, "three samples"),
            ("ttestPaired", "$Sheet1.$A$1:$A$21", {}, "one column"),
            ("anovaOneWay", "$Sheet1.$A$1:$A$21", {}, "one group"),
            ("anovaRM", "$Sheet1.$A$1:$A$21", {}, "one condition"),
            ("ttestOne", "$Sheet1.$A$1:$A$21", {"Alternative": "both"}, "unknown alternative")):
        try:
            run(command, rng, **settings)
            check(False, f"{command} refuses {what}")
        except Exception as exc:
            check(bool(str(exc)), f"{command} refuses {what}")

    sheets_before = doc.Sheets.getCount()
    for command in ("ttestOne", "ttestWelch", "ttestPaired", "anovaOneWay", "anovaWelch",
                    "anovaRM"):
        handler.trigger(command)
    check(doc.Sheets.getCount() == sheets_before, "Comparison dialogs open and cancel cleanly")
    doc.close(True)


def check_posthoc(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    rnd = random.Random(29)
    cols = [[round(rnd.gauss(10 + 1.5 * j, 1 + j), 2) for _ in range(12)] for j in range(4)]
    rows = [("Ctrl", "Low", "Mid", "High", "Group", "Value")]
    for i in range(12):
        g = ("x", "y", "z")[i % 3]
        rows.append(tuple(c[i] for c in cols) + (g, float(round(rnd.gauss(5 + i % 3, 1)))))
    rows[3] = ("n/a",) + rows[3][1:]   # A4: text, ignored
    data.getCellRangeByName("A1:F13").setDataArray(tuple(rows))
    groups = [cols[0][:2] + cols[0][3:]] + cols[1:]
    stacked = {}
    for r in rows[1:]:
        stacked.setdefault(r[4], []).append(r[5])

    def run(command, rng, **settings):
        args = (NamedValue("Command", command), NamedValue("Document", doc),
                NamedValue("Range", rng))
        return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))

    def pair(rows, a, b):
        return next(r for r in rows if r and r[0] == a and r[1] == b)

    name = run("posthocTukey", "$Sheet1.$A$1:$D$13")
    rows = report_text(doc.Sheets.getByName(name))
    ref = posthoc.tukey_hsd(groups)
    c = ref.comparisons[2]   # Ctrl - High
    row = pair(rows, "Ctrl", "High")
    check(name == "Tukey HSD" and abs(row[4] - c.statistic) < 1e-12 and abs(row[5] - c.p) < 1e-12
          and abs(row[6] - c.ci[0]) < 1e-12 and len([r for r in rows if r and r[0] in
                                                      ("Ctrl", "Low", "Mid") and r[1] in
                                                      ("Low", "Mid", "High")]) == 6,
          f"Tukey HSD: 6 pairs, Ctrl-High q = {c.statistic:.4f}, p = {c.p:.4g}")

    name = run("posthocGamesHowell", "$Sheet1.$A$1:$D$13")
    rows = report_text(doc.Sheets.getByName(name))
    c = posthoc.games_howell(groups).comparisons[5]   # Mid - High
    row = pair(rows, "Mid", "High")
    check(name == "Games-Howell" and abs(row[4] - c.statistic) < 1e-12 and
          abs(row[5] - c.df) < 1e-10 and abs(row[6] - c.p) < 1e-12,
          f"Games-Howell: Mid-High t = {c.statistic:.4f}, df = {c.df:.2f}, p = {c.p:.4g}")

    name = run("posthocDunnett", "$Sheet1.$A$1:$D$13", Control="Low", Alternative="greater")
    rows = report_text(doc.Sheets.getByName(name))
    ref = posthoc.dunnett(groups, 1, "greater")
    header = rows.index(next(r for r in rows if r[:2] == ("Group", "Difference")))
    row = find_row(rows[header:], "High")
    high = next(c for c in ref.comparisons if c.first == 3)
    check(name == "Dunnett" and find_row(rows, "Control group")[1] == "Low" and
          find_row(rows, "Ctrl") is not None and abs(row[3] - high.statistic) < 1e-12 and
          abs(row[4] - high.p) < 1e-12 and row[6] == "∞",
          f"Dunnett vs Low (greater): High t = {high.statistic:.4f}, p = {high.p:.4g}")

    name = run("posthocDunn", "$Sheet1.$F$1:$F$13", GroupRange="$Sheet1.$E$1:$E$13",
               Adjust="bonferroni")
    rows = report_text(doc.Sheets.getByName(name))
    vals = list(stacked.values())
    h, df, p = posthoc.kruskal_wallis(vals)
    c = posthoc.dunn(vals, "bonferroni").comparisons[1]   # x - z
    row = pair(rows, "x", "z")
    check(name == "Dunn" and abs(rows[rows.index(next(r for r in rows if r[0] == "H")) + 1][0] - h)
          < 1e-12 and abs(row[4] - c.statistic) < 1e-12 and abs(row[6] - c.p) < 1e-12,
          f"Dunn (Bonferroni) per group: H = {h:.4f}, x-z z = {c.statistic:.4f}, "
          f"p adj = {c.p:.4g}")

    for command, rng, settings, what in (
            ("posthocTukey", "$Sheet1.$A$1:$A$13", {}, "one group"),
            ("posthocDunnett", "$Sheet1.$A$1:$D$13", {"Control": "None"}, "an unknown control"),
            ("posthocDunn", "$Sheet1.$A$1:$D$13", {"Adjust": "hochberg"}, "an unknown adjustment")):
        try:
            run(command, rng, **settings)
            check(False, f"{command} refuses {what}")
        except Exception as exc:
            check(bool(str(exc)), f"{command} refuses {what}")

    sheets_before = doc.Sheets.getCount()
    for command in ("posthocTukey", "posthocGamesHowell", "posthocDunnett", "posthocDunn"):
        handler.trigger(command)
    check(doc.Sheets.getCount() == sheets_before, "Post hoc dialogs open and cancel cleanly")
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


def check_surface(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    # 4 x 3 grid, (1, 10) measured twice, plus a text row
    rows = [("Temp", "Time", "Yield")] + [(float(x), float(y), float(x * y))
                                         for x in (1, 2, 3, 4) for y in (10, 20, 30)]
    rows += [(1.0, 10.0, 20.0), ("n/a", 5.0, 5.0)]
    data.getCellRangeByName("A1:C15").setDataArray(tuple(rows))
    a = 2 ** 0.5
    ccd = [(-1, -1), (1, -1), (-1, 1), (1, 1), (-a, 0), (a, 0), (0, -a), (0, a), (0, 0)]
    data.getCellRangeByName("E1:G9").setDataArray(tuple((x, y, 10 - x * x - y * y)
                                                        for x, y in ccd))

    def run(rng, **settings):
        args = (NamedValue("Command", "surface"), NamedValue("Document", doc),
                NamedValue("Range", rng))
        return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))

    name = run("$Sheet1.$A$1:$C$15")
    sheet = doc.Sheets.getByName(name)
    cells = report_text(sheet)
    header = next(i for i, r in enumerate(cells) if r[0] == "Temp")
    check(name == "3D Surface" and cells[header] == ("Temp", "Time", "Yield (mean)", "n"),
          f"Surface sheet created ({name!r}) with point table")
    check(cells[header + 1] == (1.0, 10.0, 15.0, 2.0) and len(
          [r for r in cells[header + 1:] if r[0] != ""]) == 12,
          "Surface: replicates averaged, 12 points")
    check(find_row(cells, "Points plotted")[1] == 12.0 and
          find_row(cells, "Surface")[1].startswith("grid: 6 cells") and
          find_row(cells, "Rows skipped (not three numbers)")[1] == 1.0,
          "Surface: grid mesh and skipped row reported")
    page = sheet.DrawPage
    shape = page.getByIndex(0) if page.getCount() else None
    check(page.getCount() == 1 and shape.Name == "3D surface" and shape.Graphic is not None and
          shape.Graphic.MimeType == "image/svg+xml",
          "Surface: SVG picture inserted")
    check(shape.Size.Width == 16000 and shape.Size.Height == 12800, "Surface: picture size")
    check(doc.getUndoManager().getCurrentUndoActionTitle() == "LibreStats: 3D surface",
          "Surface is a single undo step")

    name = run("$Sheet1.$E$1:$G$9", HasLabels=False, Azimuth=30, Elevation=60,
               MeshLines=False, Title="CCD")
    cells = report_text(doc.Sheets.getByName(name))
    check(name == "3D Surface 2" and
          find_row(cells, "Surface")[1].startswith("Delaunay triangulation: 8 triangles") and
          find_row(cells, "View")[1] == "rotation 30°, elevation 60°" and
          next(i for i, r in enumerate(cells) if r[0] == "x") > 0,
          "Surface: scattered design triangulated, no header row")

    for rng, msg in (("$Sheet1.$A$1:$B$15", "exactly three columns"),
                     ("$Sheet1.$A$1:$C$3", "at least 3")):
        try:
            run(rng)
            check(False, f"Surface refuses {msg}")
        except Exception as exc:
            check(msg in str(exc), f"Surface refuses {msg}")
    try:
        run("$Sheet1.$A$1:$C$15", Elevation=120)
        check(False, "Surface refuses elevation 120")
    except Exception as exc:
        check("between 0 and 90" in str(exc), "Surface refuses elevation 120")

    sheets_before = doc.Sheets.getCount()
    handler.trigger("surface")
    check(doc.Sheets.getCount() == sheets_before, "Surface dialog opens and cancels cleanly")
    doc.close(True)


def summary_part(cells):
    """The summary to the right of a design table (starts with 'Design' in row 1)."""
    col = list(cells[0]).index("Design")
    return [r[col:] for r in cells]


def check_doe(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    rows = (("Name", "Unit", "Low", "High"), ("Temperature", "°C", 40.0, 80.0),
            ("Time", "min", 10.0, 30.0), ("pH", "", 5.0, 7.0), ("Stirring", "rpm", 200.0, 600.0),
            ("Catalyst", "", "Pd", "Pt"))
    data.getCellRangeByName("A1:D6").setDataArray(rows)

    def run(rng, **settings):
        args = (NamedValue("Command", "doeDesign"), NamedValue("Document", doc),
                NamedValue("Range", rng))
        return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))

    name = run("$Sheet1.$A$1:$D$5", Design="full", CentrePoints=2, Seed=11,
               Responses="Yield, Purity")
    sheet = doc.Sheets.getByName(name)
    cells = report_text(sheet)
    header = cells[0][:14]
    check(name == "Design" and header == (
        "Std order", "Run order", "Point type", "Block", "Temperature (coded)", "Time (coded)",
        "pH (coded)", "Stirring (coded)", "Temperature (°C)", "Time (min)", "pH",
        "Stirring (rpm)", "Yield", "Purity"), f"Design sheet created ({name!r}) with columns")
    runs = [r[:14] for r in cells[1:19]]
    check(len([r for r in cells[1:] if r[0] != ""]) == 18 and
          [r[1] for r in runs] == [float(i) for i in range(1, 19)] and
          sorted(r[0] for r in runs) == [float(i) for i in range(1, 19)],
          "Design: 16 + 2 centre runs in run order, standard order kept")
    check(all(r[8] == {-1.0: 40.0, 0.0: 60.0, 1.0: 80.0}[r[4]] and
              r[11] == {-1.0: 200.0, 0.0: 400.0, 1.0: 600.0}[r[7]] for r in runs),
          "Design: actual levels match coded levels")
    check(sum(1 for r in runs if r[2] == 0.0) == 2 and
          all(r[4:8] == (0.0,) * 4 for r in runs if r[2] == 0.0), "Design: centre points")
    check(all(r[12] == "" and r[13] == "" for r in runs) and
          sheet.getCellByPosition(12, 1).CellBackColor == 0xFFF5CE,
          "Design: empty, shaded response columns")
    check(find_row([r[15:] for r in cells], "Run order")[1] ==
          "randomised within blocks, seed 11", "Design: seed reported")
    again = report_text(doc.Sheets.getByName(run("$Sheet1.$A$1:$D$5", Design="full",
                                                 CentrePoints=2, Seed=11,
                                                 Responses="Yield, Purity")))
    check([r[0] for r in again[1:19]] == [r[0] for r in runs], "Design: same seed, same order")

    name = run("$Sheet1.$A$1:$D$6", Design="fractional", Randomize=False, Blocks=2)
    cells = report_text(doc.Sheets.getByName(name))
    info = [r[16:] for r in cells]
    flat = [str(v) for row in info for v in row]
    check(cells[0][14] == "Response" and cells[0][15] == "" and len([r for r in cells[1:] if r[0] != ""]) == 16 and
          find_row(info, "Design")[1] == "Fractional factorial 2^(5-1), resolution V",
          "Fractional 2^(5-1) with text factor")
    check("E = ABCD" in flat and "I = ABCDE" in flat and "Confounded with blocks" in flat,
          "Fractional: generator, defining relation and block confounding listed")
    check({r[13] for r in cells[1:17]} == {"Pd", "Pt"} and
          [r[3] for r in cells[1:17]] == sorted(r[3] for r in cells[1:17]),
          "Fractional: text levels, blocks in order")

    name = run("$Sheet1.$A$1:$D$6", Design="pb", Runs=12, Replicates=2, Blocks=2, Seed=3)
    cells = report_text(doc.Sheets.getByName(name))
    check(len([r for r in cells[1:] if r[0] != ""]) == 24 and
          any(str(v).startswith("Note: Plackett-Burman designs with a run count") for r in cells
              for v in r), "Plackett-Burman 12 runs, 2 replicates as blocks, partial aliasing note")
    check(doc.getUndoManager().getCurrentUndoActionTitle() == "LibreStats: New design",
          "Design is a single undo step")

    for settings, msg in (({"Design": "fractional", "Runs": 64}, "For 5 factors, choose"),
                          ({"Design": "full", "CentrePoints": 1}, "text levels"),
                          ({"Design": "full", "Replicates": 2, "Blocks": 3}, "multiple")):
        try:
            run("$Sheet1.$A$1:$D$6", **settings)
            check(False, f"Design refused: {msg}")
        except Exception as exc:
            check(msg in str(exc), f"Design refused: {msg}")

    # response surface designs (numeric factors only: A1:D5)
    name = run("$Sheet1.$A$1:$D$3", Design="ccd", Seed=5)
    cells = report_text(doc.Sheets.getByName(name))
    info = summary_part(cells)
    runs = [r[:10] for r in cells[1:] if r[0] != ""]
    check(len(runs) == 13 and sorted(r[2] for r in runs) == [-1.0] * 4 + [0.0] * 5 + [1.0] * 4,
          "CCD 2 factors: 4 cube, 4 axial, 5 centre points")
    check(find_row(info, "Axial distance alpha (coded)")[1] == "1.4142 (rotatable)" and
          find_row(info, "Rotatable")[1] == "yes" and
          find_row(info, "Terms in the full quadratic model")[1] == 6.0,
          "CCD: alpha, rotatability and model terms reported")
    axial_t = sorted(r[6] for r in runs if r[2] == -1.0 and r[4] != 0.0)
    check(axial_t == [31.72, 88.28], f"CCD: axial temperatures 60 -+ 20 sqrt 2, rounded "
                                     f"to 0.01 ({axial_t})")
    head = next(i for i, r in enumerate(info) if r[0] == "Factor")
    check(info[head][:8] == ("Factor", "Name", "Unit", "-alpha", "Low (-1)", "Centre (0)",
                             "High (+1)", "+alpha") and info[head + 1][4:7] == (40.0, 60.0, 80.0),
          "CCD: factor table with -alpha .. +alpha")

    name = run("$Sheet1.$A$1:$D$5", Design="ccd", Alpha="orthogonal", Blocks=3,
               CentrePoints=2, AxialCentrePoints=2, Inscribed=True, Randomize=False)
    cells = report_text(doc.Sheets.getByName(name))
    info = summary_part(cells)
    check(len([r for r in cells[1:] if r[0] != ""]) == 30 and
          find_row(info, "Blocks orthogonal to the quadratic model")[1] == "yes" and
          find_row(info, "Axial distance alpha (coded)")[1].startswith("2.0000 (orthogonal "
                                                                      "blocking, inscribed"),
          "CCD 4 factors in 3 orthogonal blocks, inscribed")
    check(min(r[8] for r in cells[1:31]) == 40.0 and max(r[8] for r in cells[1:31]) == 80.0,
          "Inscribed CCD stays within Low-High")

    name = run("$Sheet1.$A$1:$D$4", Design="bbd", Seed=9)
    cells = report_text(doc.Sheets.getByName(name))
    runs = [r for r in cells[1:] if r[0] != ""]
    check(len(runs) == 15 and {r[2] for r in runs} == {0.0, 2.0} and
          find_row(summary_part(cells), "Rotatable")[1] == "no",
          "Box-Behnken 3 factors: 12 edge midpoints + 3 centre points")
    try:
        run("$Sheet1.$A$1:$D$6", Design="ccd")
        check(False, "CCD refuses text levels")
    except Exception as exc:
        check("text levels" in str(exc), "CCD refuses text levels")

    sheets_before = doc.Sheets.getCount()
    for command in ("doeFull", "doeFractional", "doePlackettBurman", "doeCCD", "doeBoxBehnken"):
        handler.trigger(command)
    check(doc.Sheets.getCount() == sheets_before, "Design dialogs open and cancel cleanly")
    doc.close(True)


def check_screening(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    data.getCellRangeByName("A1:D5").setDataArray((
        ("Name", "Unit", "Low", "High"), ("Temperature", "°C", 20.0, 40.0),
        ("Pressure", "bar", 1.0, 2.0), ("Formaldehyde", "%", 2.0, 4.0),
        ("Stirring", "rpm", 10.0, 30.0)))

    def run(command, rng, **settings):
        args = (NamedValue("Command", command), NamedValue("Document", doc),
                NamedValue("Range", rng))
        return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))

    # Montgomery Example 6.2, randomised: responses entered by standard order
    filtration = [45, 71, 48, 65, 68, 60, 80, 65, 43, 100, 45, 104, 75, 86, 70, 96]
    name = run("doeDesign", "$Sheet1.$A$1:$D$5", Design="full", Seed=8,
               Responses="Filtration, Other")
    design = doc.Sheets.getByName(name)
    for r in range(1, 17):
        std = int(design.getCellByPosition(0, r).getValue())
        design.getCellByPosition(12, r).setValue(filtration[std - 1])
    doc.getCurrentController().setActiveSheet(design)
    report = run("doeAnalyseScreening", f"${name}.$A$1:$N$17")
    sheet = doc.Sheets.getByName(report)
    cells = report_text(sheet)
    head = next(i for i, r in enumerate(cells) if r[0] == "Term")
    effects = {r[0]: r[2] for r in cells[head + 2:head + 12]}
    check(report == "Screening analysis" and
          abs(effects["A"] - 21.625) < 1e-9 and abs(effects["AC"] + 18.125) < 1e-9 and
          abs(effects["AD"] - 16.625) < 1e-9 and len(effects) == 10,
          "Screening analysis: Montgomery 6.2 effects from a randomised design sheet")
    check(abs(find_row(cells, "Pseudo standard error (PSE)")[1] - 2.625) < 1e-12 and
          find_row(cells, "Effects beyond SME")[1] == "A, D, AC, AD",
          "Screening analysis: Lenth PSE 2.625, effects beyond SME")
    check(sheet.Charts.getCount() == 2, "Screening analysis: Pareto chart and half-normal plot")
    hn = sheet.Charts.getByName("Half-normal plot").EmbeddedObject
    series = hn.getFirstDiagram().getCoordinateSystems()[0].getChartTypes()[0].getDataSeries()[0]
    labels = [series.getDataPointByIndex(i).CustomLabelFields for i in range(15)]
    check([f[0].getString() for f in labels if f] == ["C", "D", "AD", "AC", "A"],
          "Half-normal plot: effects beyond ME labelled")
    pareto = sheet.Charts.getByName("Pareto chart").EmbeddedObject
    cs = pareto.getFirstDiagram().getCoordinateSystems()[0]
    check(cs.SwapXAndYAxis and len(cs.getChartTypes()) == 2,
          "Pareto chart: horizontal bars with a reference line")
    check(doc.getUndoManager().getCurrentUndoActionTitle() == "LibreStats: Screening analysis",
          "Screening analysis is a single undo step")

    # fraction with centre points: curvature; second response by name, main effects only
    name = run("doeDesign", "$Sheet1.$A$1:$D$5", Design="fractional", Runs=8,
               CentrePoints=3, Randomize=False)
    design = doc.Sheets.getByName(name)
    ys = [10, 14, 12, 18, 11, 15, 13, 19, 20, 21, 20.5]
    design.getCellRangeByPosition(12, 1, 12, 11).setDataArray(tuple((v,) for v in ys))
    cells = report_text(doc.Sheets.getByName(run("doeAnalyseScreening", f"${name}.$A$1:$M$12",
                                                 Model="main", Charts=False)))
    head = next(i for i, r in enumerate(cells) if r[0] == "Term")
    check([r[0] for r in cells[head + 2:head + 6]] == ["A", "B", "C", "D"] and
          find_row(cells, "Mean of centre points")[1] == 20.5 and
          str(find_row(cells, "F, p")[0]) == "F, p" and
          any(str(r[0]).startswith("Significant curvature") for r in cells),
          "Screening analysis: main effects, curvature detected")
    check(any(str(r[0]) == "A + BCD" for r in cells), "Screening analysis: alias chains listed")

    for rng, msg in (("$Sheet1.$A$1:$D$5", "design sheet"),
                     (f"${name}.$A$1:$M$12", "no response column")):
        try:
            run("doeAnalyseScreening", rng, Response="Missing" if "no response" in msg else "")
            check(False, f"Screening analysis refuses: {msg}")
        except Exception as exc:
            check(msg in str(exc), f"Screening analysis refuses: {msg}")

    # the hand-made table from documentation.md section 9.3, on another sheet at F3
    doc.Sheets.insertNewByName("Own", doc.Sheets.getCount())
    own = doc.Sheets.getByName("Own")
    table = (("Std order", "Run order", "Point type", "Block", "Temp (coded)", "Time (coded)",
              "Temp (°C)", "Time (min)", "Yield", "Purity"),
             (1, 3, 1, 1, -1, -1, 40, 10, 61.2, 97.1), (2, 1, 1, 1, 1, -1, 80, 10, 70.5, 95.4),
             (3, 5, 1, 1, -1, 1, 40, 30, 64.8, ""), (4, 2, 1, 1, 1, 1, 80, 30, 78.0, 94.2),
             (5, 4, 0, 1, 0, 0, 60, 20, 69.9, 96.0))
    own.getCellRangeByName("F3:O8").setDataArray(
        tuple(tuple(float(v) if isinstance(v, int) else v for v in r) for r in table))
    cells = report_text(doc.Sheets.getByName(run("doeAnalyseScreening", "$Own.$F$3:$O$8")))
    check(find_row(cells, "Runs analysed")[1] == 5.0 and
          [r[0] for r in cells if r[0] in ("A", "B", "AB")] == ["A", "B", "AB"] and
          find_row(cells, "Mean of centre points")[1] == 69.9,
          "Screening analysis: documented hand-made table (Yield)")
    cells = report_text(doc.Sheets.getByName(run("doeAnalyseScreening", "$Own.$F$3:$O$8",
                                                 Response="Purity", Model="main",
                                                 Charts=False)))
    check(find_row(cells, "Runs analysed")[1] == 4.0 and
          find_row(cells, "Runs without a response")[1] == 1.0 and
          [r[0] for r in cells if r[0] in ("A", "B", "AB")] == ["A", "B"],
          "Screening analysis: documented table, empty response left out (Purity)")
    own.getCellRangeByName("J8").setValue(0.4)
    try:
        run("doeAnalyseScreening", "$Own.$F$3:$O$8", Charts=False)
        check(False, "Screening analysis refuses coded value 0.4")
    except Exception as exc:
        check("must be -1, 0 or +1" in str(exc), "Screening analysis refuses coded value 0.4")

    # a CCD sheet with its response filled in (the user's report): a clear refusal, both
    # from the macro command and from the menu (dialog not opened, message instead)
    name = run("doeDesign", "$Sheet1.$A$1:$D$3", Design="ccd", Seed=1)
    ccd = doc.Sheets.getByName(name)
    n = len([r for r in report_text(ccd)[1:] if r[0] != ""])
    ccd.getCellRangeByPosition(8, 1, 8, n).setDataArray(tuple((float(i),) for i in range(n)))
    try:
        run("doeAnalyseScreening", f"${name}.$A$1:$I${n + 1}")
        check(False, "Screening analysis of a CCD: clear refusal")
    except Exception as exc:
        check("Use DOE / RSM > Analyse Response Surface Design" in str(exc),
              "Screening analysis of a CCD: clear refusal")
    doc.getCurrentController().setActiveSheet(ccd)
    sheets_before = doc.Sheets.getCount()
    handler.trigger("doeAnalyseScreening")
    check(doc.Sheets.getCount() == sheets_before, "Menu on a CCD sheet: message, no report")

    sheets_before = doc.Sheets.getCount()
    handler.trigger("doeAnalyseScreening")
    check(doc.Sheets.getCount() == sheets_before, "Screening analysis dialog opens and cancels")
    doc.close(True)


def check_rsm(desktop, handler):
    doc = desktop.loadComponentFromURL("private:factory/scalc", "_blank", 0,
                                       (PropertyValue(Name="Hidden", Value=True),))
    data = doc.Sheets.getByIndex(0)
    data.getCellRangeByName("A1:D4").setDataArray((
        ("Name", "Unit", "Low", "High"), ("Time", "min", 80.0, 90.0),
        ("Temperature", "°F", 170.0, 180.0), ("pH", "", 5.0, 7.0)))

    def run(command, rng, **settings):
        args = (NamedValue("Command", command), NamedValue("Document", doc),
                NamedValue("Range", rng))
        return handler.execute(args + tuple(NamedValue(k, v) for k, v in settings.items()))

    # Montgomery Example 11.2 typed into a randomised CCD Design sheet
    name = run("doeDesign", "$Sheet1.$A$1:$D$3", Design="ccd", Seed=3, Responses="Yield")
    design = doc.Sheets.getByName(name)
    a = 2 ** 0.5
    ys = {(-1, -1): [76.5], (-1, 1): [77.0], (1, -1): [78.0], (1, 1): [79.5],
          (0, 0): [79.9, 80.3, 80.0, 79.7, 79.8], (1, 0): [78.4], (-1, 0): [75.6],
          (0, 1): [78.5], (0, -1): [77.0]}
    for r in range(1, 14):
        key = tuple(round(design.getCellByPosition(c, r).getValue() /
                          (a if abs(abs(design.getCellByPosition(c, r).getValue()) - a) < 1e-6
                           else 1)) for c in (4, 5))
        design.getCellByPosition(8, r).setValue(ys[key].pop(0))
    doc.getCurrentController().setActiveSheet(design)
    report = run("doeAnalyseRSM", f"${name}.$A$1:$I$14")
    sheet = doc.Sheets.getByName(report)
    cells = report_text(sheet)
    head = next(i for i, r in enumerate(cells) if r[0] == "Term")
    coef = {r[0]: r[2] for r in cells[head + 1:head + 7]}
    check(report == "RSM analysis" and abs(coef["Constant"] - 79.94) < 0.001 and
          abs(coef["A^2"] + 1.376) < 0.001 and abs(coef["AB"] - 0.25) < 1e-9,
          "RSM analysis: Montgomery 11.2 coefficients from a randomised CCD sheet")
    st = next(i for i, r in enumerate(cells) if r[0] == "Stationary point and canonical analysis")
    check(abs(cells[st + 2][3] - 86.95) < 0.005 and abs(cells[st + 3][3] - 176.53) < 0.005 and
          find_row(cells, "Nature")[1] == "maximum",
          "RSM analysis: stationary point 86.95 min, 176.53 F, a maximum")
    eq = next(i for i, r in enumerate(cells) if str(r[0]).startswith("Equation in actual"))
    check(abs(cells[eq + 2][2] + 1430.52) < 0.005 and
          find_row(cells, "B")[1] == "Temperature",
          "RSM analysis: equation in actual units (constant -1430.52)")
    levels = next(i for i, r in enumerate(cells) if str(r[0]).startswith("Factor levels"))
    check(cells[levels + 2][3] == 5.0 and cells[levels + 3][3] == 5.0,
          "RSM analysis: exact half-ranges recovered from the rounded settings")
    map_cells = sheet.getCellRangeByPosition(11, 2, 31, 22)
    check(map_cells.getCellByPosition(0, 0).CellBackColor != -1 and
          sheet.getDrawPage().getCount() == 1 and
          sheet.getDrawPage().getByIndex(0).Name == "Fitted surface",
          "RSM analysis: coloured response map and 3D surface picture")
    check(doc.getUndoManager().getCurrentUndoActionTitle() ==
          "LibreStats: Response surface analysis", "RSM analysis is a single undo step")

    # Box-Behnken in 3 blocks, map of Temperature x pH
    name = run("doeDesign", "$Sheet1.$A$1:$D$4", Design="bbd", Randomize=False)
    design = doc.Sheets.getByName(name)
    n = len([r for r in report_text(design)[1:] if r[0] != ""])
    rows = report_text(design)[1:n + 1]
    design.getCellRangeByPosition(10, 1, 10, n).setDataArray(tuple(
        (60 - (r[4] - 0.3) ** 2 - 2 * r[5] ** 2 + r[6] + 0.1 * (i % 3),) for i, r in enumerate(rows)))
    cells = report_text(doc.Sheets.getByName(run("doeAnalyseRSM", f"${name}.$A$1:$K${n + 1}",
                                                 XFactor="Temperature", YFactor="C")))
    check(find_row(cells, "Runs analysed")[1] == 15.0 and
          str(cells[0][10]).startswith("Response map: predicted Response (other factors at "
                                       "their centre: Time = 85)"),
          "RSM analysis: Box-Behnken, map of chosen factors")

    # a two-level factorial with centre points: the quadratic model is refused clearly
    name = run("doeDesign", "$Sheet1.$A$1:$D$3", Design="full", CentrePoints=3,
               Randomize=False)
    design = doc.Sheets.getByName(name)
    design.getCellRangeByPosition(8, 1, 8, 7).setDataArray(tuple((float(v),) for v in
                                                                 (5, 7, 6, 9, 8, 8.2, 7.9)))
    try:
        run("doeAnalyseRSM", f"${name}.$A$1:$I$8")
        check(False, "RSM analysis refuses squares from a two-level design")
    except Exception as exc:
        check("cannot separate the square terms" in str(exc),
              "RSM analysis refuses squares from a two-level design")
    cells = report_text(doc.Sheets.getByName(run("doeAnalyseRSM", f"${name}.$A$1:$I$8",
                                                 Model="interaction", Charts=False)))
    check(find_row(cells, "Model")[1].startswith("linear + two-factor interactions"),
          "RSM analysis: interaction model on a factorial with centre points")

    sheets_before = doc.Sheets.getCount()
    handler.trigger("doeAnalyseRSM")
    check(doc.Sheets.getCount() == sheets_before, "RSM analysis dialog opens and cancels")
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
        check_surface(desktop, handler)
        check_doe(desktop, handler)
        check_screening(desktop, handler)
        check_rsm(desktop, handler)
        check_assumptions(desktop, handler)
        check_compare(desktop, handler)
        check_posthoc(desktop, handler)
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
