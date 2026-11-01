"""3-D surface plots of z over (x, y), drawn as SVG.

Calc's chart engine has no surface chart type, so the surface is drawn as a
picture. Each row of the input is one measurement (x, y, z); replicates at the
same (x, y) are averaged. Points that lie on a rectangular grid (a factorial
layout, possibly with gaps) are joined into quadrilateral tiles; scattered
points (e.g. a central composite design) are joined into triangles by a
Delaunay triangulation (Bowyer, 1981; Watson, 1981), so the surface passes
through the data and is linear between them, nothing is smoothed.

x, y and z are scaled to [-1, 1], turned by the azimuth around the vertical
axis, tilted by the elevation and projected orthographically. Tiles are drawn
from far to near (painter's algorithm, Newell et al., 1972) and coloured by
their mean height on the viridis scale (van der Walt & Smith, 2015).
"""

import math
from dataclasses import dataclass, field
from fractions import Fraction
from xml.sax.saxutils import escape

MAX_TRIANGULATED = 2000   # Bowyer-Watson here is O(n^2)
MAX_TILES = 20000
PIECES = 4000             # tiles are split into about this many pieces in total
MAX_SPLIT = 10

# viridis at 0, 1/8, ..., 1
VIRIDIS = ("#440154", "#472c7a", "#3b518b", "#2c718e", "#21908d", "#27ad81", "#5cc863",
           "#aadc32", "#fde725")

WIDTH, HEIGHT = 800, 640     # SVG size in px
FONT = "Liberation Sans, Arial, Helvetica, sans-serif"


class SurfaceError(ValueError):
    pass


def _clean(x):
    """Remove floating-point noise so 0.1 + 0.2 and 0.3 are the same grid value."""
    return float(f"{x:.12g}")


def fmt(x):
    return f"{_clean(x):.10g}" if x != 0 else "0"


@dataclass
class Surface:
    points: list                 # (x, y, mean z, number of replicates), sorted by x, y
    tiles: list                  # tuples of point indices, 4 (grid) or 3 (triangles)
    mesh: str                    # "grid" | "triangles"
    skipped: int = 0             # input rows without three numbers
    notes: list = field(default_factory=list)

    @property
    def replicated(self):
        return sum(1 for p in self.points if p[3] > 1)


def _number(v):
    if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v):
        return float(v)
    return None


def build_surface(rows):
    """Surface from rows of (x, y, z); rows without three numbers are skipped."""
    sums, skipped = {}, 0
    for row in rows:
        values = [_number(v) for v in list(row)[:3]]
        if len(values) < 3 or None in values:
            if any(v is not None and str(v).strip() for v in row):
                skipped += 1
            continue
        key = (_clean(values[0]), _clean(values[1]))
        total, n = sums.get(key, (0.0, 0))
        sums[key] = (total + values[2], n + 1)
    if len(sums) < 3:
        raise SurfaceError("A surface needs at least 3 different (x, y) points with "
                           "numeric x, y and z.")
    points = [(x, y, total / n, n) for (x, y), (total, n) in sorted(sums.items())]
    if len({p[0] for p in points}) < 2 or len({p[1] for p in points}) < 2:
        raise SurfaceError("x and y must each take at least two different values.")

    surface = Surface(points, [], "grid", skipped)
    quads = grid_tiles(points)
    if quads is not None:
        surface.tiles = quads
    else:
        if len(points) > MAX_TRIANGULATED:
            raise SurfaceError(f"The points do not form a grid, and triangulating more than "
                               f"{MAX_TRIANGULATED} scattered points is too slow.")
        surface.mesh = "triangles"
        surface.tiles = delaunay(_scaled_xy(points))
        if not surface.tiles:
            raise SurfaceError("All (x, y) points lie on one line, so they span no surface.")
    if len(surface.tiles) > MAX_TILES:
        raise SurfaceError(f"The surface would have {len(surface.tiles)} tiles; the limit is "
                           f"{MAX_TILES}. Use a coarser grid.")
    if surface.mesh == "grid":
        xs, ys = {p[0] for p in points}, {p[1] for p in points}
        missing = (len(xs) - 1) * (len(ys) - 1) - len(surface.tiles)
        if missing:
            surface.notes.append(f"{missing} grid cell(s) lack a corner point and are left "
                                 "open.")
    if surface.replicated:
        surface.notes.append(f"{surface.replicated} (x, y) point(s) were measured more than "
                             "once; their z values are averaged.")
    if len({p[2] for p in points}) == 1:
        surface.notes.append("All z values are equal: the surface is flat.")
    return surface


def grid_tiles(points):
    """Quadrilaterals between neighbouring grid values, or None if some point
    would belong to no tile (the points are then treated as scattered)."""
    xs = sorted({p[0] for p in points})
    ys = sorted({p[1] for p in points})
    index = {(p[0], p[1]): i for i, p in enumerate(points)}
    tiles, used = [], set()
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            corners = ((xs[i], ys[j]), (xs[i + 1], ys[j]), (xs[i + 1], ys[j + 1]),
                       (xs[i], ys[j + 1]))
            if all(c in index for c in corners):
                tile = tuple(index[c] for c in corners)
                tiles.append(tile)
                used.update(tile)
    return tiles if len(used) == len(points) else None


def _scaled_xy(points):
    return [(_scale(p[0], *_extent([q[0] for q in points])),
             _scale(p[1], *_extent([q[1] for q in points]))) for p in points]


def _extent(values):
    return min(values), max(values)


def _scale(v, lo, hi):
    return 0.0 if hi == lo else 2 * (v - lo) / (hi - lo) - 1


# --- Delaunay triangulation (Bowyer-Watson) ------------------------------------

def _incircle(a, b, c, d):
    """> 0 if d lies inside the circle through a, b, c (counter-clockwise).
    Rounding errors near 0 are settled with exact rational arithmetic."""
    def det(a, b, c, d):
        adx, ady = a[0] - d[0], a[1] - d[1]
        bdx, bdy = b[0] - d[0], b[1] - d[1]
        cdx, cdy = c[0] - d[0], c[1] - d[1]
        ad, bd, cd = adx * adx + ady * ady, bdx * bdx + bdy * bdy, cdx * cdx + cdy * cdy
        terms = (ad * (bdx * cdy - cdx * bdy), bd * (cdx * ady - adx * cdy),
                 cd * (adx * bdy - bdx * ady))
        return sum(terms), sum(abs(t) for t in terms)

    value, size = det(a, b, c, d)
    if abs(value) > 1e-9 * size:
        return value
    exact = [tuple(Fraction(v) for v in p) for p in (a, b, c, d)]
    return det(*exact)[0]


def delaunay(xy):
    """Counter-clockwise triangles (index triples) of the Delaunay triangulation."""
    big = 1e6
    pts = list(xy) + [(-big, -big), (big, -big), (0.0, big)]
    n = len(xy)
    triangles = {(n, n + 1, n + 2)}
    for i in range(n):
        p = pts[i]
        bad = [t for t in triangles if _incircle(pts[t[0]], pts[t[1]], pts[t[2]], p) > 0]
        edges = {(t[k], t[(k + 1) % 3]) for t in bad for k in range(3)}
        triangles.difference_update(bad)
        # cavity boundary: edges whose twin is not in the cavity; they keep their orientation
        triangles.update((a, b, i) for a, b in edges if (b, a) not in edges)
    return sorted(t for t in triangles if max(t) < n)


# --- drawing --------------------------------------------------------------------

def nice_ticks(lo, hi, target=5):
    """Round tick values (1, 2 or 5 times a power of ten) between lo and hi;
    the step whose number of ticks is closest to target wins."""
    if hi == lo:
        return [lo]
    mag = 10 ** math.floor(math.log10((hi - lo) / target))
    best = None
    for step in (m * mag * f for f in (0.1, 1, 10) for m in (1, 2, 5)):
        first = math.ceil(lo / step - 1e-9)
        count = math.floor(hi / step + 1e-9) - first + 1
        if best is None or abs(count - target) <= abs(best[0] - target):
            best = (count, step, first)
    count, step, first = best
    return [_clean((first + k) * step) for k in range(count)]


def color(t):
    """Viridis colour for t in [0, 1]."""
    t = min(max(t, 0.0), 1.0) * (len(VIRIDIS) - 1)
    i = min(int(t), len(VIRIDIS) - 2)
    f = t - i
    a, b = VIRIDIS[i], VIRIDIS[i + 1]
    rgb = [round(int(a[k:k + 2], 16) * (1 - f) + int(b[k:k + 2], 16) * f) for k in (1, 3, 5)]
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def _darker(hex_color, factor=0.7):
    return "#{:02x}{:02x}{:02x}".format(*(round(int(hex_color[k:k + 2], 16) * factor)
                                          for k in (1, 3, 5)))


class View:
    """Orthographic camera: turn by azimuth a about the z axis, tilt by elevation e."""

    def __init__(self, azimuth, elevation):
        a, e = math.radians(azimuth), math.radians(elevation)
        self.ca, self.sa, self.ce, self.se = math.cos(a), math.sin(a), math.cos(e), math.sin(e)

    def project(self, x, y, z):
        """Screen position (u right, v up) and distance d from the camera."""
        xr = x * self.ca - y * self.sa
        yr = x * self.sa + y * self.ca
        return xr, z * self.ce + yr * self.se, yr * self.ce - z * self.se

    def depth(self, x, y, z):
        return self.project(x, y, z)[2]


def _mix(points, weights):
    return tuple(sum(w * p[k] for p, w in zip(points, weights)) for k in range(len(points[0])))


def subdivide(corners, m):
    """Split a tile into m x m pieces lying on the same surface: linear
    interpolation within a triangle, bilinear within a grid quadrilateral.
    Yields (piece corners, edges of the piece on the tile's outline)."""
    if len(corners) == 3:
        a, b, c = corners

        def at(i, j):
            return _mix((a, b, c), ((m - i - j) / m, i / m, j / m))

        for i in range(m):
            for j in range(m - i):
                edges = [e for e, on in (((0, 1), j == 0), ((1, 2), i + j + 1 == m),
                                         ((2, 0), i == 0)) if on]
                yield [at(i, j), at(i + 1, j), at(i, j + 1)], edges
                if i + j + 2 <= m:
                    yield [at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)], []
        return

    def at(i, j):
        s, t = i / m, j / m
        return _mix(corners, ((1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t))

    for i in range(m):
        for j in range(m):
            edges = [e for e, on in (((0, 1), j == 0), ((1, 2), i == m - 1),
                                     ((2, 3), j == m - 1), ((3, 0), i == 0)) if on]
            yield [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)], edges


def _text(x, y, s, size=12, anchor="middle", bold=False, rotate=None, fill="#333333"):
    extra = ' font-weight="bold"' if bold else ""
    if rotate is not None:
        extra += f' transform="rotate({rotate} {x:.1f} {y:.1f})"'
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
            f'text-anchor="{anchor}" fill="{fill}"{extra}>{escape(str(s))}</text>')


def _line(p, q, stroke="#999999", width=1.0):
    return (f'<line x1="{p[0]:.1f}" y1="{p[1]:.1f}" x2="{q[0]:.1f}" y2="{q[1]:.1f}" '
            f'stroke="{stroke}" stroke-width="{width}"/>')


def _polygon(corners, fill, stroke, width=0.6):
    pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in corners)
    return (f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{width}" '
            'stroke-linejoin="round"/>')


def _text_width(s, size):
    return 0.55 * size * len(str(s))


def render_svg(surface, azimuth=-60.0, elevation=30.0, title="", x_title="x", y_title="y",
               z_title="z", mesh_lines=True):
    """The surface as an SVG document (str)."""
    if not 0 <= elevation <= 90:
        raise SurfaceError("The elevation must be between 0 and 90 degrees.")
    view = View(azimuth, elevation)
    pts = surface.points
    ranges = [_extent([p[k] for p in pts]) for k in range(3)]
    ticks = [nice_ticks(*r) for r in ranges]
    scaled = [tuple(_scale(p[k], *ranges[k]) for k in range(3)) for p in pts]

    def sc(k, v):
        return _scale(v, *ranges[k])

    # front edges for the x and y tick labels, back panes, z axis on the leftmost edge
    yf = -1.0 if view.depth(0, -1, -1) < view.depth(0, 1, -1) else 1.0
    xf = -1.0 if view.depth(-1, 0, -1) < view.depth(1, 0, -1) else 1.0
    zc = min(((cx, cy) for cx in (-1, 1) for cy in (-1, 1)),
             key=lambda c: (round(view.project(c[0], c[1], 0)[0], 9), -view.depth(c[0], c[1], 0)))

    out = 1.12   # tick labels this far outside the box
    title_out = 1.45
    anchors = [(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    anchors += [(0, yf * title_out, -1), (xf * title_out, 0, -1)]
    proj = [view.project(*a) for a in anchors]
    umin, umax = min(p[0] for p in proj), max(p[0] for p in proj)
    vmin, vmax = min(p[1] for p in proj), max(p[1] for p in proj)

    left, right, top, bottom = 90, 150, 60, 50
    s = min((WIDTH - left - right) / (umax - umin), (HEIGHT - top - bottom) / (vmax - vmin))
    ox = left + ((WIDTH - left - right) - s * (umax - umin)) / 2
    oy = top + ((HEIGHT - top - bottom) - s * (vmax - vmin)) / 2

    def screen(x, y, z):
        u, v, _ = view.project(x, y, z)
        return ox + s * (u - umin), oy + s * (vmax - v)

    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
           f'viewBox="0 0 {WIDTH} {HEIGHT}">',
           f'<rect x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="#ffffff"/>']
    if title:
        svg.append(_text(WIDTH / 2, 32, title, size=17, bold=True, fill="#222222"))

    # back panes with grid lines
    xb, yb = -xf, -yf
    panes = [
        ([(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1)],
         [((sc(0, t), -1, -1), (sc(0, t), 1, -1)) for t in ticks[0]] +
         [((-1, sc(1, t), -1), (1, sc(1, t), -1)) for t in ticks[1]]),
        ([(xb, -1, -1), (xb, 1, -1), (xb, 1, 1), (xb, -1, 1)],
         [((xb, sc(1, t), -1), (xb, sc(1, t), 1)) for t in ticks[1]] +
         [((xb, -1, sc(2, t)), (xb, 1, sc(2, t))) for t in ticks[2]]),
        ([(-1, yb, -1), (1, yb, -1), (1, yb, 1), (-1, yb, 1)],
         [((sc(0, t), yb, -1), (sc(0, t), yb, 1)) for t in ticks[0]] +
         [((-1, yb, sc(2, t)), (1, yb, sc(2, t))) for t in ticks[2]]),
    ]
    for corners, grid in panes:
        svg.append(_polygon([screen(*c) for c in corners], "#f4f4f4", "#c8c8c8", 1.0))
        svg += [_line(screen(*a), screen(*b), "#dddddd", 0.8) for a, b in grid]

    # tiles, split into pieces for smooth colours, far to near
    zlo, zhi = ranges[2]
    m = max(1, min(MAX_SPLIT, math.isqrt(PIECES // len(surface.tiles))))
    pieces = []
    for tile in surface.tiles:
        corners = [scaled[i] + (pts[i][2],) for i in tile]
        for piece, edges in subdivide(corners, m):
            depth = sum(view.depth(*c[:3]) for c in piece) / len(piece)
            mean_z = sum(c[3] for c in piece) / len(piece)
            pieces.append((depth, [screen(*c[:3]) for c in piece], edges,
                           color(0.5 if zhi == zlo else (mean_z - zlo) / (zhi - zlo))))
    pieces.sort(key=lambda p: p[0], reverse=True)
    for _, corners, edges, fill in pieces:
        svg.append(_polygon(corners, fill, fill, 0.8))   # stroke hides seams between pieces
        if mesh_lines:
            svg += [_line(corners[a], corners[b], _darker(fill), 0.7) for a, b in edges]

    # x and y axes on the front floor edges
    def axis(k, edge_point, outward, values, name):
        a = screen(*edge_point(-1))
        b = screen(*edge_point(1))
        svg.append(_line(a, b, "#555555", 1.2))
        centre = screen(0, 0, -1)
        for t in values:
            p = edge_point(sc(k, t))
            q = tuple(c + o * (out - 1) for c, o in zip(p, outward))
            p2, q2 = screen(*p), screen(*q)
            svg.append(_line(p2, q2, "#555555", 1.0))
            du = q2[0] - p2[0]
            anchor = "start" if du > 2 else "end" if du < -2 else "middle"
            svg.append(_text(q2[0], q2[1] + 12, fmt(t), size=11, anchor=anchor))
        mid = edge_point(0)
        m = screen(*(c + o * (title_out - 1) for c, o in zip(mid, outward)))
        du = m[0] - centre[0]
        anchor = "start" if du > 0.15 * WIDTH else "end" if du < -0.15 * WIDTH else "middle"
        svg.append(_text(m[0], m[1] + 16, name, size=13, anchor=anchor, bold=True))

    axis(0, lambda v: (v, yf, -1), (0, yf, 0), ticks[0], x_title)
    axis(1, lambda v: (xf, v, -1), (xf, 0, 0), ticks[1], y_title)

    # z axis on the leftmost vertical edge
    bottom_pt, top_pt = screen(zc[0], zc[1], -1), screen(zc[0], zc[1], 1)
    svg.append(_line(bottom_pt, top_pt, "#555555", 1.2))
    widest = 0
    for t in ticks[2]:
        x, y = screen(zc[0], zc[1], sc(2, t))
        svg.append(_line((x, y), (x - 6, y), "#555555", 1.0))
        svg.append(_text(x - 9, y + 4, fmt(t), size=11, anchor="end"))
        widest = max(widest, _text_width(fmt(t), 11))
    zx = min(bottom_pt[0], top_pt[0]) - 18 - widest
    svg.append(_text(zx, (bottom_pt[1] + top_pt[1]) / 2, z_title, size=13, bold=True,
                     rotate=-90))

    # colour legend
    lx, ly, lw, lh = WIDTH - 110, top + 30, 18, HEIGHT - top - bottom - 60
    steps = 64
    for i in range(steps):
        y0 = ly + lh * (1 - (i + 1) / steps)
        fill = color((i + 0.5) / steps)
        svg.append(f'<rect x="{lx}" y="{y0:.2f}" width="{lw}" height="{lh / steps + 0.6:.2f}" '
                   f'fill="{fill}" stroke="none"/>')
    svg.append(f'<rect x="{lx}" y="{ly}" width="{lw}" height="{lh}" fill="none" '
               'stroke="#888888" stroke-width="0.8"/>')
    for t in ticks[2]:
        y = ly + lh * (0.5 if zhi == zlo else 1 - (t - zlo) / (zhi - zlo))
        svg.append(_line((lx + lw, y), (lx + lw + 5, y), "#555555", 1.0))
        svg.append(_text(lx + lw + 8, y + 4, fmt(t), size=11, anchor="start"))
    svg.append(_text(lx + lw / 2, ly - 12, z_title, size=12, bold=True))

    svg.append("</svg>")
    return "\n".join(svg)
