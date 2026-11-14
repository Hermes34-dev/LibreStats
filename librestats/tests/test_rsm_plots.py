"""Check librestats.rsm_plots: the plotted values (traces are model predictions, they meet
at the anchor, boundary optima sit at the line ends), the SVG content (points, labels,
ticks, column means) and the greyscale style; and the greyscale/footnote options of the
3D surface."""

import os
import re
import sys
import unittest
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import doe, rsm, rsm_plots, surface  # noqa: E402

NS = "{http://www.w3.org/2000/svg}"
Y = [24.28989323, 33.66123277, 39.83253039, 18.65134664, 12.52643918, 24.52400157,
     19.43370304, 18.04370825, 11.58758592, 19.25931338, 16.25878132, 9.714146726,
     13.70971019, 33.93868664, 22.13181921, 51.13647479, 26.76090827, 37.48920235,
     15.37816958, 31.05661974, 11.95990485, 18.31395824, 17.97089328, 11.73768574,
     8.318061785, 12.1663204, 14.31427388, 10.81211525]


def design2():
    d = doe.build_rsm([doe.Factor(n, "", 0.0, 1.0) for n in "TtE"], "cccd", randomize=False)
    x = [r.coded for r in d.runs]
    return x, rsm.fit_model(x, [1] * 28, Y)


def parse(svg):
    root = ET.fromstring(svg)
    texts = ["".join(t.itertext()) for t in root.iter(NS + "text")]
    return root, texts


class TestTrace(unittest.TestCase):
    def test_traces_are_predictions_and_meet_at_the_anchor(self):
        x, fit = design2()
        a = fit.maximum[0]
        svg, traces, y_anchor = rsm_plots.optimum_trace(
            fit, a, ["Time (min)", "Temp (C)", "EtOH (%)"], [13, 125, 50],
            [7.135, 44.6, 29.73], "Optimum trace: Yield", "full quadratic: A, B", "Yield")
        self.assertAlmostEqual(y_anchor, fit.maximum[1], places=9)
        for i, (xs, ys) in enumerate(traces):
            lo, hi = fit.bounds[i]
            self.assertAlmostEqual(xs[0], lo)
            self.assertAlmostEqual(xs[-1], hi)
            for xv, yv in zip(xs[::10], ys[::10]):
                z = list(a)
                z[i] = xv
                self.assertAlmostEqual(yv, fit.predict(z), places=12)
            # the trace passes through the anchor (all traces meet there)
            z = list(a)
            self.assertAlmostEqual(max(ys), y_anchor, places=6)   # maximum of each trace
        root, texts = parse(svg)
        # one open circle per factor at the anchor height (plus 3 in the legend)
        circles = [c for c in root.iter(NS + "circle") if c.get("fill") == "#ffffff"]
        self.assertEqual(len(circles), 6)
        cy = {round(float(c.get("cy")), 2) for c in circles[:3]}
        self.assertEqual(len(cy), 1)
        # a boundary optimum: the circle sits at the end of its line
        lines = list(root.iter(NS + "polyline"))
        for i, c in enumerate(circles[:3]):
            pts = [tuple(map(float, p.split(","))) for p in lines[i].get("points").split()]
            if abs(a[i] - fit.bounds[i][1]) < 1e-9:
                self.assertAlmostEqual(float(c.get("cx")), pts[-1][0], places=1)
            if abs(a[i] - fit.bounds[i][0]) < 1e-9:
                self.assertAlmostEqual(float(c.get("cx")), pts[0][0], places=1)
        # stacked real-unit axes with their own ticks
        for t in ("Time (min)", "Temp (C)", "EtOH (%)", "100", "200", "25"):
            self.assertIn(t, texts)
        self.assertIn("Optimum trace: Yield", texts)

    def test_centre_anchor_is_a_perturbation_plot(self):
        x, fit = design2()
        svg, traces, y0 = rsm_plots.optimum_trace(fit, [0.0] * 3, ["A", "B", "C"], [0] * 3,
                                                  [1] * 3, "Perturbation plot", "", "y",
                                                  "centre")
        self.assertAlmostEqual(y0, fit.predict([0, 0, 0]))
        self.assertIn("circles: centre, predicted", " ".join(parse(svg)[1]))


class TestObservedPredicted(unittest.TestCase):
    def test_points_labels_and_square_scale(self):
        x, fit = design2()
        obs = [d.observed for d in fit.diagnostics]
        pred = [d.fitted for d in fit.diagnostics]
        labels = [str(i + 1) if d.outlier or d.influential else None
                  for i, d in enumerate(fit.diagnostics)]
        svg = rsm_plots.observed_vs_predicted(obs, pred, labels, fit.r2, fit.r2_pred,
                                              "Observed vs predicted Yield", "m", "Yield")
        root, texts = parse(svg)
        self.assertEqual(len(list(root.iter(NS + "circle"))), 28)
        self.assertIn("3", texts)                                 # run 3: the outlier
        self.assertIn(f"R² = {fit.r2:.3f};  predicted R² = {fit.r2_pred:.3f}", texts)
        frame = [r for r in root.iter(NS + "rect") if r.get("stroke") == "#888888"][0]
        self.assertEqual(frame.get("width"), frame.get("height"))
        diag = [ln for ln in root.iter(NS + "line") if ln.get("stroke-dasharray") == "8,5"][0]
        # y = x runs corner to corner of the square frame
        self.assertAlmostEqual(float(diag.get("x2")) - float(diag.get("x1")),
                               float(diag.get("y1")) - float(diag.get("y2")), places=2)


class TestResidualsByFactor(unittest.TestCase):
    def test_columns_ticks_means(self):
        x, fit = design2()
        segments = []
        for j, name in enumerate(("Time (min)", "Temp (C)", "EtOH (%)")):
            levels = sorted({round(v[j], 9) for v in x})
            segments.append((name, [(lv, f"{lv:.3g}") for lv in levels],
                             [levels.index(round(v[j], 9)) for v in x]))
        res = [d.residual for d in fit.diagnostics]
        svg = rsm_plots.residuals_by_factor(segments, res, "Residuals", "m", "Residual")
        root, texts = parse(svg)
        self.assertEqual(len(list(root.iter(NS + "circle"))), 3 * 28)
        self.assertEqual(sum(len(s[1]) for s in segments), 15)          # 5 levels each
        means = [ln for ln in root.iter(NS + "line") if ln.get("stroke-width") == "2.4"]
        self.assertEqual(len(means), 15)
        for name in ("Time (min)", "Temp (C)", "EtOH (%)"):
            self.assertIn(name, texts)
        self.assertEqual(svg, rsm_plots.residuals_by_factor(segments, res, "Residuals", "m",
                                                            "Residual"))   # same jitter
        # each mean bar sits at its column's mean residual: map y back to residual units
        # with the y-axis tick labels (text right-aligned left of the axis)
        ticks = [(float(t.get("y")) - 4, float("".join(t.itertext())))
                 for t in root.iter(NS + "text")
                 if t.get("text-anchor") == "end" and re.fullmatch(r"-?[\d.]+",
                                                                   "".join(t.itertext()))]
        (y1, v1), (y2, v2) = ticks[0], ticks[-1]
        value = lambda y: v1 + (y - y1) * (v2 - v1) / (y2 - y1)  # noqa: E731
        expected = []
        for _, levels, index in segments:
            for li in range(len(levels)):
                col = [r for r, idx in zip(res, index) if idx == li]
                expected.append(sum(col) / len(col))
        for bar, mean in zip(means, expected):
            self.assertAlmostEqual(value(float(bar.get("y1"))), mean, delta=0.01)

    def test_greyscale_has_no_colours(self):
        x, fit = design2()
        segments = [("A", [(-1, "lo"), (1, "hi")], [0 if v[0] < 0 else 1 for v in x])]
        svg = rsm_plots.residuals_by_factor(segments, [d.residual for d in fit.diagnostics],
                                            "t", "s", "r", greyscale=True)
        for colour in rsm_plots.PALETTE:
            self.assertNotIn(colour, svg)


class TestSurfaceOptions(unittest.TestCase):
    def test_greyscale_and_footnote(self):
        rows = [(a, b, a * a + b) for a in range(-2, 3) for b in range(-2, 3)]
        s = surface.build_surface(rows)
        svg = surface.render_svg(s, title="Fitted y", greyscale=True, subtitle="model",
                                 footnote="C fixed at the predicted maximum = 100")
        self.assertIn("C fixed at the predicted maximum = 100", svg)
        self.assertIn(">model<", svg)
        fills = set(re.findall(r'fill="(#[0-9a-f]{6})"', svg))
        self.assertFalse(fills & set(surface.VIRIDIS))
        self.assertEqual(surface.grey(0), "#2e2e2e")
        self.assertEqual(surface.grey(1), "#ebebeb")       # 255 x 0.92 = 234.6

    def test_pretty(self):
        self.assertEqual(rsm_plots.pretty("A^2B, C^3"), "A²B, C³")


if __name__ == "__main__":
    unittest.main()
