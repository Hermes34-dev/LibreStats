"""Check librestats.surface: mesh building against scipy's Delaunay, projection and SVG."""

import math
import os
import sys
import unittest
import xml.etree.ElementTree as ET

import numpy as np
from scipy.spatial import Delaunay

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "pythonpath"))
from librestats import surface as sf  # noqa: E402

RNG = np.random.default_rng(11)
SVG = "{http://www.w3.org/2000/svg}"


def triangle_set(triangles):
    return {frozenset(t) for t in triangles}


def area(xy, t):
    (ax, ay), (bx, by), (cx, cy) = (xy[i] for i in t)
    return ((bx - ax) * (cy - ay) - (cx - ax) * (by - ay)) / 2


class TestDelaunay(unittest.TestCase):
    def test_matches_scipy_random_points(self):
        for n in (3, 10, 60, 300):
            xy = [tuple(p) for p in RNG.uniform(-1, 1, (n, 2))]
            ref = triangle_set(Delaunay(np.array(xy)).simplices.tolist())
            self.assertEqual(triangle_set(sf.delaunay(xy)), ref, f"n = {n}")

    def test_counter_clockwise_and_covers_hull(self):
        xy = [tuple(p) for p in RNG.normal(0, 1, (80, 2))]
        tris = sf.delaunay(xy)
        self.assertTrue(all(area(xy, t) > 0 for t in tris))
        hull = Delaunay(np.array(xy)).convex_hull   # total area = hull area
        ref = sum(abs(area(xy, t)) for t in Delaunay(np.array(xy)).simplices.tolist())
        self.assertAlmostEqual(sum(area(xy, t) for t in tris), ref, delta=1e-9)
        self.assertGreater(len(hull), 2)

    def test_cocircular_points(self):
        # rotatable central composite design: the 8 outer points lie on one circle
        a = math.sqrt(2)
        xy = [(-1, -1), (1, -1), (-1, 1), (1, 1), (-a, 0), (a, 0), (0, -a), (0, a), (0, 0)]
        tris = sf.delaunay(xy)
        self.assertEqual(len(tris), 8)              # fan around the centre
        self.assertTrue(all(8 in t for t in tris))
        self.assertAlmostEqual(sum(area(xy, t) for t in tris), 4 * a, delta=1e-12)
        # square grid: every cell is cocircular, still a valid triangulation
        grid = [(float(i), float(j)) for i in range(5) for j in range(5)]
        tris = sf.delaunay(grid)
        self.assertEqual(len(tris), 32)
        self.assertAlmostEqual(sum(area(grid, t) for t in tris), 16, delta=1e-12)

    def test_collinear_gives_nothing(self):
        self.assertEqual(sf.delaunay([(0, 0), (1, 1), (2, 2), (3, 3)]), [])


class TestBuild(unittest.TestCase):
    def test_full_grid_uses_quads(self):
        rows = [(x, y, x * y) for x in (1, 2, 3, 4) for y in (10, 20, 30)]
        s = sf.build_surface(rows)
        self.assertEqual((s.mesh, len(s.tiles), len(s.points)), ("grid", 6, 12))
        self.assertEqual(s.notes, [])

    def test_grid_with_gap_and_replicates(self):
        rows = [(x, y, x + y) for x in (0, 1, 2) for y in (0, 1, 2) if (x, y) != (2, 2)]
        rows += [(0, 0, 4.0), ("a", 1, 2), (None, None, None), (1, "", 3)]
        s = sf.build_surface(rows)
        self.assertEqual(s.mesh, "grid")
        self.assertEqual(len(s.tiles), 3)
        self.assertEqual(s.points[0], (0.0, 0.0, 2.0, 2))      # (0 + 4) / 2
        self.assertEqual(s.skipped, 2)                          # the all-empty row is not counted
        self.assertTrue(any("left open" in n for n in s.notes))
        self.assertTrue(any("averaged" in n for n in s.notes))

    def test_float_noise_still_a_grid(self):
        # 0.1 * 3 = 0.30000000000000004 must meet 0.3 typed in another row
        rows = [(0.1 * i, y, i * y) for i in range(4) for y in (0, 1)]
        rows[-1] = (0.3, 1, 3)
        s = sf.build_surface(rows)
        self.assertEqual((s.mesh, len(s.tiles)), ("grid", 3))

    def test_scattered_points_are_triangulated(self):
        a = math.sqrt(2)
        ccd = [(-1, -1), (1, -1), (-1, 1), (1, 1), (-a, 0), (a, 0), (0, -a), (0, a), (0, 0)]
        s = sf.build_surface([(x, y, x * x + y) for x, y in ccd])
        self.assertEqual((s.mesh, len(s.tiles)), ("triangles", 8))

    def test_errors(self):
        with self.assertRaises(sf.SurfaceError):
            sf.build_surface([(1, 1, 1), (2, 2, 2)])
        with self.assertRaises(sf.SurfaceError):
            sf.build_surface([(1, y, y) for y in range(5)])               # one x value
        with self.assertRaises(sf.SurfaceError):
            sf.build_surface([(i, i, i) for i in range(5)])                 # on a line


class TestDrawing(unittest.TestCase):
    def test_nice_ticks(self):
        self.assertEqual(sf.nice_ticks(0, 10), [0, 2, 4, 6, 8, 10])
        self.assertEqual(sf.nice_ticks(60.02, 84.97), [65, 70, 75, 80])
        self.assertEqual(sf.nice_ticks(41.3, 98.7), [50, 60, 70, 80, 90])
        self.assertEqual(sf.nice_ticks(0.1, 0.35), [0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
        self.assertEqual(sf.nice_ticks(-1, 1), [-1, -0.5, 0, 0.5, 1])
        self.assertEqual(sf.nice_ticks(5, 5), [5])

    def test_viridis_ends(self):
        self.assertEqual(sf.color(0), "#440154")
        self.assertEqual(sf.color(1), "#fde725")
        self.assertEqual(sf.color(0.5), "#21908d")

    def test_projection(self):
        v = sf.View(0, 0)                      # looking along +y, level
        self.assertEqual(v.project(1, 0, 0), (1, 0, 0))
        u, h, d = v.project(0, 1, 0)
        self.assertAlmostEqual(d, 1)           # further away
        v = sf.View(90, 0)                     # turned a quarter: +x points away
        self.assertAlmostEqual(v.project(1, 0, 0)[2], 1)
        v = sf.View(0, 90)                     # from straight above
        u, h, d = v.project(0, 1, 1)
        self.assertAlmostEqual(h, 1)           # +y is up on the screen
        self.assertAlmostEqual(d, -1)          # high points are nearer
        v = sf.View(-37, 23)                   # rotation keeps lengths
        for p in ((1, 0, 0), (0.3, -0.7, 0.5)):
            self.assertAlmostEqual(sum(c * c for c in v.project(*p)), sum(c * c for c in p))

    def test_subdivide_stays_on_surface(self):
        tri = [(0, 0, 0, 1.0), (1, 0, 0, 2.0), (0, 1, 1, 3.0)]
        pieces = list(sf.subdivide(tri, 4))
        self.assertEqual(len(pieces), 16)
        total = sum(abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2
                    for (a, b, c), _ in pieces)
        self.assertAlmostEqual(total, 0.5)
        for piece, _ in pieces:
            for x, y, z, w in piece:            # on the plane z = y, w = 1 + x + 2y
                self.assertAlmostEqual(z, y)
                self.assertAlmostEqual(w, 1 + x + 2 * y)
        self.assertEqual(sum(len(e) for _, e in pieces), 12)   # 3 sides x 4 segments

        quad = [(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 1, 1), (0, 1, 0, 0)]  # z = x * y
        pieces = list(sf.subdivide(quad, 3))
        self.assertEqual(len(pieces), 9)
        for piece, _ in pieces:
            for x, y, z, _w in piece:
                self.assertAlmostEqual(z, x * y)
        self.assertEqual(sum(len(e) for _, e in pieces), 12)

    def svg(self, azimuth, elevation, **kw):
        rows = [(x, y, math.sin(x) * math.cos(y)) for x in np.linspace(0, 3, 8)
                for y in np.linspace(0, 3, 6)]
        return ET.fromstring(sf.render_svg(sf.build_surface(rows), azimuth, elevation,
                                           "T & <t>", **kw))

    def test_svg_is_valid_for_all_angles(self):
        for az in range(-180, 181, 30):
            for el in (0, 15, 45, 90):
                root = self.svg(az, el)
                self.assertEqual(root.tag, SVG + "svg")
                coords = [float(v) for p in root.iter(SVG + "polygon")
                          for xy in p.get("points").split() for v in xy.split(",")]
                self.assertTrue(all(-1 <= c <= 801 for c in coords), (az, el))
        texts = [t.text for t in root.iter(SVG + "text")]
        self.assertIn("T & <t>", texts)
        self.assertIn("x", texts)

    def test_far_pieces_drawn_first(self):
        # a flat plane seen from the front: pieces must go from the back (top) to the front
        rows = [(x, y, 0.0) for x in range(3) for y in range(3)]
        root = ET.fromstring(sf.render_svg(sf.build_surface(rows), 0, 30, mesh_lines=False))
        tiles = [p for p in root.iter(SVG + "polygon") if p.get("fill") == sf.color(0.5)]
        ys = [sum(float(xy.split(",")[1]) for xy in p.get("points").split()) for p in tiles]
        self.assertEqual(ys, sorted(ys))

    def test_bad_elevation(self):
        with self.assertRaises(sf.SurfaceError):
            self.svg(0, 95)


if __name__ == "__main__":
    unittest.main()
