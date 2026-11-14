"""Diagnostic and interpretation plots of a fitted response surface model, drawn as SVG
(the Calc chart engine cannot draw stacked axes or per-segment category axes).

- Observed vs predicted: one point per run, the line y = x, a square plot with the
  same scale on both axes, R^2 and predicted R^2; flagged runs labelled.
- Residuals by factor: one shared residual axis, the x axis split into one segment
  per factor, one column of points per design level (low to high, labelled with
  the real setting), a small deterministic jitter and the mean of each column. A
  curve across a factor's columns suggests a missing term in that factor; spread
  that changes across levels suggests unequal variance.
- Optimum trace: for each factor, the predicted response as that factor moves over
  its range with the others held at the anchor (the optimum, or the centre for the
  classic perturbation plot; Myers et al., 2016). One shared position scale, with a
  real-unit axis per factor stacked under the plot; the anchor is marked on every
  trace, so the traces meet at the anchor's predicted response.

Colour plots use the LibreOffice chart palette in a fixed factor order; the
greyscale (journal) style uses black and grey with line styles and marker shapes.
"""

import math
import random
from xml.sax.saxutils import escape

from .surface import FONT, fmt, nice_ticks

# LibreOffice's default chart colours, without the yellow (unreadable on white)
PALETTE = ("#004586", "#ff420e", "#579d1c", "#7e0021", "#83caff", "#314004", "#aecf00")
GREY_PALETTE = ("#000000", "#555555", "#000000", "#555555", "#000000", "#555555", "#000000")
DASHES = (None, "9,5", "2,4", "10,4,2,4", "16,5", "5,5", "2,2")   # per factor
ACCENT = "#c9211e"        # reference lines, flagged runs (colour style)
TICK = "#555555"
GRID = "#e4e4e4"


def pretty(name):
    """'A^2B' -> 'A²B' for plot text."""
    return name.replace("^2", "²").replace("^3", "³")


def palette(greyscale):
    return GREY_PALETTE if greyscale else PALETTE


class Svg:
    def __init__(self, width, height):
        self.width, self.height = width, height
        self.items = [f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>']

    def text(self, x, y, s, size=12, anchor="middle", bold=False, rotate=None,
             fill="#333333", italic=False):
        extra = (' font-weight="bold"' if bold else "") + \
            (' font-style="italic"' if italic else "")
        rot = f' transform="rotate({rotate} {x:.1f} {y:.1f})"' if rotate else ""
        self.items.append(
            f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
            f'text-anchor="{anchor}" fill="{fill}"{extra}{rot}>{escape(str(s))}</text>')

    def line(self, p, q, stroke=TICK, width=1.0, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.items.append(f'<line x1="{p[0]:.2f}" y1="{p[1]:.2f}" x2="{q[0]:.2f}" '
                          f'y2="{q[1]:.2f}" stroke="{stroke}" stroke-width="{width}"{d}/>')

    def polyline(self, points, stroke, width=2.0, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        pts = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        self.items.append(f'<polyline points="{pts}" fill="none" stroke="{stroke}" '
                          f'stroke-width="{width}" stroke-linejoin="round"{d}/>')

    def circle(self, x, y, r, fill, stroke=None, width=1.0):
        s = f' stroke="{stroke}" stroke-width="{width}"' if stroke else ""
        self.items.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}"{s}/>')

    def rect(self, x, y, w, h, fill="none", stroke=TICK, width=1.0):
        self.items.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
                          f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')

    def render(self):
        return "\n".join([f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" '
                          f'height="{self.height}" viewBox="0 0 {self.width} '
                          f'{self.height}">'] + self.items + ["</svg>"])


def _wrap(text, width):
    """Split text into lines of at most width characters, at ', ' or spaces."""
    words, lines, line = text.split(" "), [], ""
    for w in words:
        if line and len(line) + 1 + len(w) > width:
            lines.append(line)
            line = w
        else:
            line = f"{line} {w}" if line else w
    return lines + ([line] if line else [])


def _titles(svg, title, subtitle, width, chars=95):
    """Title and wrapped subtitle; returns the y below them."""
    svg.text(width / 2, 28, title, size=17, bold=True, fill="#222222")
    y = 28
    for line in _wrap(pretty(subtitle), chars) if subtitle else []:
        y += 17
        svg.text(width / 2, y, line, size=12, fill="#444444")
    return y + 14


def _range(values, pad=0.05):
    lo, hi = min(values), max(values)
    if hi == lo:
        lo, hi = lo - 1, hi + 1
    span = hi - lo
    return lo - pad * span, hi + pad * span


def _y_axis(svg, x0, top, height, lo, hi, title, ticks=None):
    ticks = ticks or [t for t in nice_ticks(lo, hi) if lo - 1e-12 <= t <= hi + 1e-12]

    def ys(v):
        return top + height * (hi - v) / (hi - lo)
    svg.line((x0, top), (x0, top + height), TICK, 1.2)
    widest = 0
    for t in ticks:
        y = ys(t)
        svg.line((x0 - 5, y), (x0, y), TICK, 1.0)
        svg.text(x0 - 8, y + 4, fmt(t), size=11, anchor="end")
        widest = max(widest, 0.55 * 11 * len(fmt(t)))
    svg.text(x0 - 16 - widest, top + height / 2, title, size=13, bold=True, rotate=-90)
    return ys, ticks


# --- 1. observed vs predicted --------------------------------------------------------

def observed_vs_predicted(observed, predicted, labels, r2, r2_pred, title, subtitle,
                          response, greyscale=False):
    """Square plot of observed (y) against predicted (x). labels: text per run for the
    runs to label (outliers, influential), None for the others."""
    W = 640
    svg = Svg(W, 10)
    top = _titles(svg, title, subtitle, W, chars=80) + 10
    left, size = 95, W - 95 - 40
    svg.height = top + size + 58
    svg.items[0] = f'<rect x="0" y="0" width="{W}" height="{svg.height}" fill="#ffffff"/>'
    lo, hi = _range(list(observed) + list(predicted))
    ticks = [t for t in nice_ticks(lo, hi) if lo <= t <= hi]

    def sx(v):
        return left + size * (v - lo) / (hi - lo)

    def sy(v):
        return top + size * (hi - v) / (hi - lo)
    for t in ticks:                                         # grid
        svg.line((sx(t), top), (sx(t), top + size), GRID, 0.8)
        svg.line((left, sy(t)), (left + size, sy(t)), GRID, 0.8)
    svg.rect(left, top, size, size, stroke="#888888", width=1.0)
    svg.line((sx(lo), sy(lo)), (sx(hi), sy(hi)), "#000000" if greyscale else ACCENT, 1.4,
             "8,5")
    for t in ticks:
        svg.line((sx(t), top + size), (sx(t), top + size + 5), TICK)
        svg.text(sx(t), top + size + 18, fmt(t), size=11)
        svg.line((left - 5, sy(t)), (left, sy(t)), TICK)
        svg.text(left - 8, sy(t) + 4, fmt(t), size=11, anchor="end")
    svg.text(left + size / 2, top + size + 40, f"Predicted {response}", size=13, bold=True)
    svg.text(left - 62, top + size / 2, f"Observed {response}", size=13, bold=True,
             rotate=-90)
    colour = palette(greyscale)[0]
    for o, p, lab in zip(observed, predicted, labels):
        if lab:
            svg.circle(sx(p), sy(o), 5.5, "#ffffff", "#000000" if greyscale else ACCENT, 1.8)
            svg.text(sx(p) + 8, sy(o) - 6, lab, size=12, anchor="start", bold=True)
        else:
            svg.circle(sx(p), sy(o), 4.5, colour, "#ffffff", 0.6)
    r2p = "n/a" if r2_pred is None else f"{r2_pred:.3f}"
    svg.text(left + 10, top + 20, f"R² = {r2:.3f};  predicted R² = {r2p}",
             size=13, anchor="start", fill="#222222")
    svg.text(left + size - 8, top + size - 10, "dashed: observed = predicted", size=11,
             anchor="end", fill="#666666", italic=True)
    if any(labels):
        svg.text(left + size - 8, top + size - 26, "open circles: outlier or influential "
                 "(std order)", size=11, anchor="end", fill="#666666", italic=True)
    return svg.render()


# --- 2. residuals by factor ----------------------------------------------------------

def residuals_by_factor(segments, residuals, title, subtitle, y_title, greyscale=False,
                        guides=False, seed=7):
    """segments: [(factor label, [(level value, tick text)], [level index per run])];
    residuals: one value per run (None = not available). guides: dotted lines at +-2."""
    k = len(segments)
    W = max(760, 60 * sum(len(s[1]) for s in segments) + 160)
    H = 560
    svg = Svg(W, H)
    top = _titles(svg, title, subtitle, W) + 12
    left, right, bottom = 90, 25, 95
    height = H - top - bottom
    values = [r for r in residuals if r is not None]
    m = max([abs(v) for v in values] + ([2.2] if guides else [1e-12]))
    lo, hi = -1.1 * m, 1.1 * m
    ys, _ = _y_axis(svg, left, top, height, lo, hi, y_title)
    width = (W - left - right) / k
    colours = palette(greyscale)
    svg.line((left, ys(0)), (W - right, ys(0)), "#000000", 1.0, "6,4")
    if guides:
        for g in (-2, 2):
            svg.line((left, ys(g)), (W - right, ys(g)), "#888888", 0.9, "2,3")
    rng = random.Random(seed)
    for j, (label, levels, index) in enumerate(segments):
        x0 = left + j * width
        if j:
            svg.line((x0, top), (x0, top + height), "#bbbbbb", 1.0)
        n = len(levels)
        step = width / (n + 1)
        jitter = min(0.28 * step, 12)
        colour = colours[j % len(colours)]
        for li, (_, tick) in enumerate(levels):
            cx = x0 + step * (li + 1)
            svg.line((cx, top + height), (cx, top + height + 5), TICK)
            svg.text(cx, top + height + 18, tick, size=11)
            col = [r for r, idx in zip(residuals, index) if idx == li and r is not None]
            for r in col:
                svg.circle(cx + rng.uniform(-jitter, jitter), ys(r), 4.0, colour,
                           "#ffffff", 0.6)
            if col:
                mean = sum(col) / len(col)
                svg.line((cx - 1.6 * jitter - 4, ys(mean)), (cx + 1.6 * jitter + 4, ys(mean)),
                         "#000000", 2.4)
        svg.text(x0 + width / 2, top + height + 42, label, size=13, bold=True)
    svg.line((left, top + height), (W - right, top + height), TICK, 1.2)
    svg.text(W - right, H - 14, "short bars: mean of each column; points jittered "
             "sideways", size=11, anchor="end", fill="#666666", italic=True)
    return svg.render()


# --- 4. optimum trace ----------------------------------------------------------------

def optimum_trace(fit, anchor, factors, centres, halves, title, subtitle, response,
                  anchor_name="optimum", greyscale=False, points=61):
    """Trace of each factor over its coded range (fit.bounds) with the others at the
    anchor (coded point). factors: axis labels (with units)."""
    k = fit.k
    W, plot_h, axis_h = 820, 380, 46
    traces = []
    for i in range(k):
        lo, hi = fit.bounds[i]
        xs = [lo + (hi - lo) * j / (points - 1) for j in range(points)]
        ys = []
        for x in xs:
            z = list(anchor)
            z[i] = x
            ys.append(fit.predict(z))
        traces.append((xs, ys))
    y_anchor = fit.predict(list(anchor))
    allv = [v for _, ys in traces for v in ys]
    lo_y, hi_y = _range(allv, 0.08)

    legend_rows = math.ceil(k / 3)
    svg = Svg(W, 10)              # height set below
    top = _titles(svg, title, subtitle, W) + 10
    left, right = 95, 30
    pw = W - left - right
    H = top + plot_h + 20 + k * axis_h + 26 + legend_rows * 22 + 16
    svg.height = H
    svg.items[0] = f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>'
    ys_map, ticks = _y_axis(svg, left, top, plot_h, lo_y, hi_y, f"Predicted {response}")
    for t in ticks:
        svg.line((left, ys_map(t)), (left + pw, ys_map(t)), GRID, 0.8)
    svg.rect(left, top, pw, plot_h, stroke="#888888", width=1.0)
    svg.line((left, ys_map(y_anchor)), (left + pw, ys_map(y_anchor)), "#aaaaaa", 0.9, "3,4")
    colours = palette(greyscale)

    def ux(i, x):
        lo, hi = fit.bounds[i]
        return left + pw * (0.5 if hi == lo else (x - lo) / (hi - lo))
    for i, (xs, ys) in enumerate(traces):
        svg.polyline([(ux(i, x), ys_map(y)) for x, y in zip(xs, ys)],
                     colours[i % len(colours)], 2.4, DASHES[i % len(DASHES)])
    for i in range(k):                       # anchor marks on top of all lines
        svg.circle(ux(i, anchor[i]), ys_map(y_anchor), 6.5, "#ffffff",
                   colours[i % len(colours)], 2.2)
    # stacked real-unit axes, one per factor
    y = top + plot_h + 20
    for i in range(k):
        lo, hi = fit.bounds[i]
        real_lo, real_hi = sorted((centres[i] + halves[i] * lo, centres[i] + halves[i] * hi))
        colour = colours[i % len(colours)]
        svg.line((left, y), (left + pw, y), colour, 1.6, DASHES[i % len(DASHES)])
        for t in nice_ticks(real_lo, real_hi):
            if not real_lo - 1e-9 <= t <= real_hi + 1e-9:
                continue
            coded = (t - centres[i]) / halves[i]
            x = ux(i, coded)
            svg.line((x, y), (x, y + 5), TICK)
            svg.text(x, y + 17, fmt(t), size=11)
        svg.text(left - 8, y + 5, factors[i], size=12, anchor="end", bold=True, fill=colour
                 if not greyscale else "#222222")
        y += axis_h
    # legend: line style per factor
    y += 6
    per_row = 3
    col_w = pw / per_row
    for i in range(k):
        r, c = divmod(i, per_row)
        x = left + c * col_w
        yy = y + r * 22
        colour = colours[i % len(colours)]
        svg.line((x, yy), (x + 40, yy), colour, 2.4, DASHES[i % len(DASHES)])
        svg.circle(x + 20, yy, 5, "#ffffff", colour, 2.0)
        svg.text(x + 48, yy + 4, factors[i], size=12, anchor="start")
    svg.text(left + pw / 2, top - 6, f"circles: {anchor_name}, predicted "
             f"{fmt(round(y_anchor, 4))} (all traces meet there)", size=11, fill="#666666",
             italic=True)
    return svg.render(), traces, y_anchor
