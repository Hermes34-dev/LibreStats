"""Two-level screening designs: full and fractional factorials, Plackett-Burman.

Factors are coded -1 (low) and +1 (high); centre points are coded 0. Runs are
listed in standard (Yates) order, with factor A changing fastest.

- Full factorial 2^k: all 2^k combinations.
- Fractional factorial 2^(k-p): the minimum-aberration fraction (Fries &
  Hunter, 1980) from ff_table; the defining relation, resolution (Box &
  Hunter, 1961) and alias chains are reported.
- Plackett-Burman (1946): N = 8, 12, 16, 20 or 24 runs from the cyclic
  generator rows of the paper, plus a row of all minus signs.

Blocks split each replicate by confounding interactions with blocks. The block
generators are chosen to confound no main effect and as few low-order
interactions as possible (lexicographically smallest pattern of the orders of
the confounded effects, aliases included). Run order is randomised within
blocks.
"""

import math
import random
from dataclasses import dataclass, field
from collections import Counter
from itertools import combinations, combinations_with_replacement

from .ff_table import DESIGNS as FF_DESIGNS

LETTERS = "ABCDEFGHJKLMNOPQRSTUVWXYZ"   # I is the identity in alias lists
KINDS = ("full", "fractional", "pb")
MAX_FULL = 8
MAX_FRACTIONAL = 11
MAX_RUNS = 5000
ALIAS_ORDER = 3          # interactions above this order are left out of alias lists

# Plackett & Burman (1946): first row of each cyclic design
PB_ROWS = {
    8: "+++-+--",
    12: "++-+++---+-",
    16: "++++-+-++--+---",
    20: "++--++++-+-+----++-",
    24: "+++++-+-++--++--+-+----",
}


class DesignError(ValueError):
    pass


def _clean(x):
    return float(f"{x:.12g}")


@dataclass
class Factor:
    name: str
    unit: str
    low: object              # number, or text for a categorical factor
    high: object

    @property
    def numeric(self):
        return not isinstance(self.low, str)

    def actual(self, coded, scale=1.0):
        """Actual level of a coded value; scale = alpha puts +-alpha at the low and high
        levels (inscribed central composite design)."""
        if not self.numeric:
            return self.low if coded < 0 else self.high
        centre, half = (self.low + self.high) / 2, (self.high - self.low) / 2
        return _clean(centre + coded * half / scale)

    def setting(self, coded, scale=1.0):
        """actual() rounded to a practical precision: about 4 significant digits of the
        half-range (e.g. 948.8 rpm for 200-600 rpm)."""
        value = self.actual(coded, scale)
        if not self.numeric:
            return value
        half = abs(self.high - self.low) / 2
        return _clean(round(value, max(0, 3 - math.floor(math.log10(half)))))

    def label(self):
        return f"{self.name} ({self.unit})" if self.unit else self.name


def _cell(v):
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return float(v) if math.isfinite(v) else None
    v = str(v).strip()
    return v or None


def read_factors(rows):
    """Factors from rows of (name, unit, low, high) or (name, low, high)."""
    factors = []
    for row in rows:
        cells = [_cell(v) for v in row]
        if all(c is None for c in cells):
            continue
        if len(cells) == 4:
            name, unit, low, high = cells
        elif len(cells) == 3:
            (name, low, high), unit = cells, None
        else:
            raise DesignError("The factor table needs the columns Name, Unit, Low, High "
                              "(or Name, Low, High).")
        if name is None:
            raise DesignError("Every factor needs a name.")
        name = name if isinstance(name, str) else f"{name:g}"
        if low is None or high is None:
            raise DesignError(f"Factor {name}: enter both a low and a high level.")
        if isinstance(low, str) != isinstance(high, str):
            raise DesignError(f"Factor {name}: the levels must both be numbers or both be "
                              "text.")
        if low == high:
            raise DesignError(f"Factor {name}: the low and high levels are equal.")
        unit = "" if unit is None else unit if isinstance(unit, str) else f"{unit:g}"
        factors.append(Factor(name, unit, low, high))
    names = [f.name for f in factors]
    if len(set(names)) < len(names):
        raise DesignError("Factor names must be different.")
    if len(factors) < 2:
        raise DesignError("A screening design needs at least 2 factors.")
    if len(factors) > len(LETTERS):
        raise DesignError(f"At most {len(LETTERS)} factors.")
    return factors


# --- effects as bit masks (bit i = factor i) --------------------------------------

def effect_name(mask):
    return "".join(LETTERS[i] for i in range(len(LETTERS)) if mask >> i & 1) or "I"


def _order(mask):
    return bin(mask).count("1")


def _sorted_effects(masks):
    return sorted(masks, key=lambda m: (_order(m), [i for i in range(32) if m >> i & 1]))


def _mask(word):
    return sum(1 << LETTERS.index(c) for c in word)


def standard_order(n):
    """2^n runs in Yates order (first factor changes fastest)."""
    return [[1 if r >> i & 1 else -1 for i in range(n)] for r in range(1 << n)]


def subgroup(generators):
    """All non-identity products of the generator words."""
    words = []
    for w in generators:
        words += [w] + [w ^ v for v in words]
    return words


@dataclass
class Fraction:
    k: int
    runs: int
    resolution: int = None           # None for a full factorial
    generators: list = field(default_factory=list)   # (added factor, product mask)
    words: list = field(default_factory=list)        # defining relation without I

    def matrix(self):
        n = self.runs.bit_length() - 1
        rows = standard_order(n)
        for factor, product in self.generators:
            for row in rows:
                row.append(math.prod(row[i] for i in range(n) if product >> i & 1))
        return rows

    def aliases(self, mask):
        return {mask ^ w for w in self.words} | {mask}


def full_factorial(k):
    if not 2 <= k <= MAX_FULL:
        raise DesignError(f"A full factorial design here has 2 to {MAX_FULL} factors "
                          f"({1 << MAX_FULL} runs).")
    return Fraction(k, 1 << k)


def fractional_choices(k):
    """(runs, resolution) of the available fractions for k factors."""
    return sorted((runs, res) for (kk, runs), (res, _, _) in FF_DESIGNS.items() if kk == k)


def roman(n):
    return ("", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X")[n]


def describe_choices(choices):
    return ", ".join(f"{runs} runs (resolution {roman(res)})" for runs, res in choices)


def fractional_factorial(k, runs=None):
    """Minimum-aberration 2^(k-p) fraction; runs=None: the smallest of resolution IV
    or more if there is one, else the largest fraction."""
    choices = fractional_choices(k)
    if not choices:
        raise DesignError(f"Fractional factorial designs here have 3 to {MAX_FRACTIONAL} "
                          "factors.")
    if runs is None:
        runs = next((r for r, res in choices if res >= 4), choices[-1][0])
    if runs not in dict(choices):
        raise DesignError(f"For {k} factors, choose {describe_choices(choices)}.")
    res, gens, _ = FF_DESIGNS[(k, runs)]
    n = runs.bit_length() - 1
    generators = [(n + j, _mask(g)) for j, g in enumerate(gens)]
    words = subgroup([product | 1 << factor for factor, product in generators])
    return Fraction(k, runs, res, generators, words)


def pb_choices(k):
    return [n for n in PB_ROWS if n >= k + 1]


def plackett_burman(k, runs=None):
    choices = pb_choices(k)
    if not choices:
        raise DesignError(f"Plackett-Burman designs here have up to {max(PB_ROWS) - 1} "
                          "factors.")
    runs = choices[0] if runs is None else runs
    if runs not in choices:
        raise DesignError(f"For {k} factors, choose {', '.join(map(str, choices))} runs.")
    first = [1 if c == "+" else -1 for c in PB_ROWS[runs]]
    m = runs - 1
    rows = [[first[(j - i) % m] for j in range(m)] for i in range(m)] + [[-1] * m]
    return [row[:k] for row in rows]


# --- blocking -----------------------------------------------------------------------

def _pattern(fraction, masks):
    """Count of confounded effects by order (aliases included): (order 1, order 2, ...)."""
    counts = [0] * (fraction.k + 1)
    for m in masks:
        for e in fraction.aliases(m):
            counts[_order(e)] += 1
    return counts[1:]


def block_generators(fraction, blocks):
    """Effects (masks of basic factors) whose signs define the blocks."""
    b = blocks.bit_length() - 1
    n = fraction.runs.bit_length() - 1
    if b == 0:
        return []
    single = {c: _pattern(fraction, [c]) for c in range(1, 1 << n)}
    cands = sorted((c for c in single if single[c][0] == 0), key=lambda c: (single[c], c))
    best = [None, None]

    def extend(chosen, elements, start):
        pattern = _pattern(fraction, elements)
        if pattern[0] or (best[0] is not None and pattern > best[0]):
            return
        if len(chosen) == b:
            if best[0] is None or pattern < best[0]:
                best[0], best[1] = pattern, list(chosen)
            return
        for i in range(start, len(cands)):
            c = cands[i]
            if c in elements:
                continue
            extend(chosen + [c], elements + [c] + [c ^ e for e in elements], i + 1)

    extend([], [], 0)
    if best[1] is None:
        raise DesignError(f"{blocks} blocks would confound a main effect with blocks. Use "
                          "fewer blocks or more runs.")
    return best[1]


# --- the design ---------------------------------------------------------------------

CUBE, CENTRE, AXIAL, EDGE = 1, 0, -1, 2      # point types (as Minitab)


@dataclass
class Run:
    std: int
    block: int
    ptype: int
    coded: list
    order: int = 0

    @property
    def centre(self):
        return self.ptype == CENTRE


@dataclass
class Design:
    kind: str
    factors: list
    base_runs: int
    runs: list                       # in run order
    replicates: int
    blocks: int
    centre_points: int               # per block
    fraction: Fraction = None        # full / fractional factorials
    block_effects: list = field(default_factory=list)   # masks confounded with blocks
    seed: int = None                 # None: not randomised
    notes: list = field(default_factory=list)
    alpha: float = None              # central composite: axial distance (coded)
    scale: float = 1.0               # alpha for an inscribed design, else 1
    axial_centre_points: int = None  # central composite in blocks: centre points, axial block
    properties: list = field(default_factory=list)   # (name, value) rows for the summary

    def title(self):
        k = len(self.factors)
        if self.kind == "full":
            return f"Full factorial 2^{k}"
        if self.kind == "fractional":
            f = self.fraction
            return f"Fractional factorial 2^({k}-{k - (f.runs.bit_length() - 1)}), " \
                   f"resolution {roman(f.resolution)}"
        if self.kind == "ccd":
            f = self.fraction
            cube = f"2^{k}" if f.runs == 1 << k else \
                f"2^({k}-{k - (f.runs.bit_length() - 1)})"
            kind = "inscribed central composite" if self.scale != 1 else "Central composite"
            return f"{kind[0].upper() + kind[1:]} ({cube} + {2 * k} axial points)"
        if self.kind == "bbd":
            return f"Box-Behnken, {k} factors"
        return f"Plackett-Burman, {self.base_runs} runs"

    def generator_text(self):
        return [f"{LETTERS[factor]} = {effect_name(product)}"
                for factor, product in self.fraction.generators]

    def defining_relation(self):
        return "I = " + " = ".join(effect_name(w) for w in _sorted_effects(self.fraction.words))

    def alias_chains(self, max_order=ALIAS_ORDER):
        """Alias chains of the main effects and two-factor interactions, each once."""
        k = len(self.factors)
        effects = [1 << i for i in range(k)] + [1 << i | 1 << j for i, j in
                                                 combinations(range(k), 2)]
        chains, seen = [], set()
        for e in effects:
            if e in seen:
                continue
            chain = self.fraction.aliases(e)
            seen |= chain
            shown = [m for m in _sorted_effects(chain) if _order(m) <= max_order]
            chains.append(" + ".join(effect_name(m) for m in shown))
        return chains

    def block_text(self):
        return [" + ".join(effect_name(m) for m in _sorted_effects(self.fraction.aliases(e))
                           if _order(m) <= ALIAS_ORDER) or effect_name(e)
                for e in _sorted_effects(subgroup(self.block_effects))]


def _run_order(runs, n_blocks, randomize, seed):
    """Shuffle the runs within each block (blocks stay in order); number them."""
    if randomize:
        if seed is None:
            seed = random.SystemRandom().randrange(1, 1000000)
        rng = random.Random(seed)
        ordered = []
        for b in range(1, n_blocks + 1):
            group = [r for r in runs if r.block == b]
            rng.shuffle(group)
            ordered += group
        runs = ordered
    else:
        seed = None
    for i, r in enumerate(runs, 1):
        r.order = i
    return runs, seed


def build_design(factors, kind, runs=None, centre_points=0, replicates=1, blocks=1,
                 randomize=True, seed=None):
    k = len(factors)
    if kind not in KINDS:
        raise DesignError(f"Unknown design {kind!r}.")
    if centre_points < 0 or replicates < 1 or blocks < 1:
        raise DesignError("Centre points must be 0 or more; replicates and blocks 1 or more.")
    if centre_points and not all(f.numeric for f in factors):
        raise DesignError("Centre points need numeric levels for every factor; "
                          + ", ".join(f.name for f in factors if not f.numeric)
                          + " has text levels.")

    fraction, notes = None, []
    if kind == "full":
        fraction = full_factorial(k)
        matrix = fraction.matrix()
    elif kind == "fractional":
        fraction = fractional_factorial(k, runs)
        matrix = fraction.matrix()
    else:
        matrix = plackett_burman(k, runs)
        n = len(matrix)
        if n & (n - 1):
            notes.append("Plackett-Burman designs with a run count that is not a power of 2 "
                         "have partial aliasing: each two-factor interaction is partly "
                         "aliased with every main effect it does not involve.")
        else:
            notes.append(f"With {n} runs this Plackett-Burman design is a regular "
                         "resolution III fraction: main effects are aliased with two-factor "
                         "interactions.")
        if n - 1 > k:
            notes.append(f"{n - 1 - k} column(s) of the design are unused; they can serve to "
                         "estimate the error.")
    base_runs = len(matrix)

    if blocks > 1 and blocks % replicates:
        raise DesignError(f"With {replicates} replicate(s), the number of blocks must be a "
                          "multiple of the replicates: each replicate is one block or is split "
                          "into 2, 4 or 8 blocks.")
    per_rep = blocks // replicates if blocks > 1 else 1
    if per_rep not in (1, 2, 4, 8):
        raise DesignError("Each replicate can be split into 1, 2, 4 or 8 blocks.")
    if per_rep > 1 and fraction is None:
        raise DesignError("Plackett-Burman designs can only be blocked by replicate (blocks = "
                          "replicates).")
    if per_rep > 1 and base_runs // per_rep < 2:
        raise DesignError(f"{base_runs} runs cannot be split into {per_rep} blocks.")
    block_effects = block_generators(fraction, per_rep) if per_rep > 1 else []

    def block_of(row):
        signs = [math.prod(row[i] for i in range(k) if e >> i & 1) for e in block_effects]
        return tuple(signs)

    # blocks within a replicate, numbered in order of first appearance in standard order
    labels = {}
    for row in matrix:
        labels.setdefault(block_of(row), len(labels))
    total = base_runs * replicates + centre_points * blocks
    if total > MAX_RUNS:
        raise DesignError(f"The design would have {total} runs; the limit is {MAX_RUNS}.")

    design_runs, std = [], 0
    n_blocks = blocks if blocks > 1 else 1
    for rep in range(replicates):
        groups = [[] for _ in range(per_rep)]
        for row in matrix:
            groups[labels[block_of(row)]].append(row)
        for j, group in enumerate(groups):
            block = rep * per_rep + j + 1 if blocks > 1 else 1
            for row in group:
                std += 1
                design_runs.append(Run(std, block, CUBE, list(row)))
            for _ in range(centre_points):
                std += 1
                design_runs.append(Run(std, block, CENTRE, [0] * k))

    design_runs, seed = _run_order(design_runs, n_blocks, randomize, seed)

    if centre_points and kind == "pb":
        notes.append("Centre points do not change the Plackett-Burman design; they give a "
                     "check for curvature and an estimate of pure error.")
    return Design(kind, factors, base_runs, design_runs, replicates, blocks, centre_points,
                  fraction, block_effects, seed, notes)


# --- response surface designs -------------------------------------------------------

RSM_KINDS = ("ccd", "bbd")
ALPHAS = ("rotatable", "orthogonal", "face", "custom")
# Box & Hunter (1957): lambda4 of uniform-precision rotatable designs, k = 2..7
UNIFORM_PRECISION = {2: 0.7844, 3: 0.8385, 4: 0.8704, 5: 0.8918, 6: 0.9070, 7: 0.9184}
CCD_CUBE_BLOCKS = {1: 1, 2: 1, 3: 2, 5: 4}     # blocks per replicate -> blocks of the cube
BLOCK_CENTRE = 2                               # default centre points per block

# Box & Behnken (1960): sets of factors run as 2^m factorials with the others at 0
_PAIRS = {k: tuple(combinations(range(k), 2)) for k in (3, 4, 5)}
BB_SETS = {**_PAIRS,
           6: ((0, 1, 3), (1, 2, 4), (2, 3, 5), (0, 3, 4), (1, 4, 5), (0, 2, 5)),
           7: ((3, 4, 5), (0, 5, 6), (1, 4, 6), (0, 1, 3), (2, 3, 6), (0, 2, 4), (1, 2, 5))}
BB_CENTRE = {3: 3, 4: 3, 5: 6, 6: 6, 7: 6}
BB_BLOCKS = {4: (((0, 1), (2, 3)), ((0, 3), (1, 2)), ((0, 2), (1, 3)))}   # orthogonal blocks


def ccd_cube(k, fraction=True):
    """Factorial part of a central composite design: the full 2^k, or for 5 or more
    factors (fraction=True) the smallest fraction of resolution V or higher."""
    if not 2 <= k <= max(UNIFORM_PRECISION):
        raise DesignError(f"Central composite designs here have 2 to {max(UNIFORM_PRECISION)} "
                          "factors.")
    if fraction and k >= 5:
        runs = next(r for r, res in fractional_choices(k) if res >= 5)
        return fractional_factorial(k, runs)
    return full_factorial(k)


def rotatable_alpha(n_cube):
    return n_cube ** 0.25


def uniform_precision_centre(k, n_cube):
    """Centre points of the uniform-precision rotatable design (Box & Hunter, 1957):
    N = lambda4 (nF + 2 alpha^2)^2 / nF with alpha^2 = sqrt(nF)."""
    n = UNIFORM_PRECISION[k] * (n_cube + 2 * math.sqrt(n_cube)) ** 2 / n_cube
    return max(1, round(n - n_cube - 2 * k))


def orthogonal_alpha(k, n_cube, cube_centre, axial_centre):
    """Axial distance that makes the axial block orthogonal to the cube block(s):
    nF / (nF + n_cf) = 2 alpha^2 / (2k + n_ca) (Box & Hunter, 1957)."""
    return math.sqrt(n_cube * (2 * k + axial_centre) / (2 * (n_cube + cube_centre)))


def model_row(x):
    """Full quadratic model: 1, linear, squares, two-factor interactions."""
    return ([1.0] + list(x) + [v * v for v in x] +
            [x[i] * x[j] for i, j in combinations(range(len(x)), 2)])


def model_terms(k):
    return 1 + 2 * k + k * (k - 1) // 2


def is_rotatable(points, tol=1e-9):
    """Moment conditions of rotatability up to order 4 (Box & Hunter, 1957): odd
    moments vanish, [ii], [iiii] and [iijj] are the same for all factors, and
    [iiii] = 3 [iijj]."""
    k, n = len(points[0]), len(points)

    def moment(idx):
        return sum(math.prod(p[i] for i in idx) for p in points) / n

    scale = max(moment((i, i, i, i)) for i in range(k)) or 1.0
    for d in range(1, 5):
        for idx in combinations_with_replacement(range(k), d):
            if any(idx.count(i) % 2 for i in set(idx)) and abs(moment(idx)) > tol * scale:
                return False
    second = [moment((i, i)) for i in range(k)]
    fourth = [moment((i, i, i, i)) for i in range(k)]
    mixed = [moment((i, i, j, j)) for i, j in combinations(range(k), 2)]
    same = lambda v: max(v) - min(v) <= tol * scale  # noqa: E731
    return same(second) and same(fourth) and same(mixed) and \
        abs(fourth[0] - 3 * mixed[0]) <= tol * scale


def blocks_orthogonal(runs, tol=1e-9):
    """True if the block effects are orthogonal to every term of the quadratic model:
    each model column has the same mean in every block. None for a single block."""
    blocks = sorted({r.block for r in runs})
    if len(blocks) < 2:
        return None
    rows = [model_row(r.coded)[1:] for r in runs]
    overall = [sum(c) / len(rows) for c in zip(*rows)]
    scale = max(map(abs, overall)) or 1.0
    for b in blocks:
        sub = [row for row, r in zip(rows, runs) if r.block == b]
        means = [sum(c) / len(sub) for c in zip(*sub)]
        if any(abs(m - o) > tol * max(scale, 1.0) for m, o in zip(means, overall)):
            return False
    return True


def box_behnken_sets(k):
    if k not in BB_SETS:
        raise DesignError(f"Box-Behnken designs here have 3 to {max(BB_SETS)} factors.")
    return BB_SETS[k]


def _subfactorial(k, factors):
    rows = []
    for signs in standard_order(len(factors)):
        row = [0] * k
        for f, s in zip(factors, signs):
            row[f] = s
        rows.append(row)
    return rows


def _per_replicate_blocks(blocks, replicates, allowed, what):
    if blocks == 1:
        return 1
    if blocks % replicates or blocks // replicates not in allowed:
        choices = ", ".join(str(a) for a in allowed)
        raise DesignError(f"{what}: with {replicates} replicate(s) the number of blocks must "
                          f"be the replicates times {choices}.")
    return blocks // replicates


def build_rsm(factors, kind, alpha="rotatable", custom_alpha=None, fraction=True,
              inscribed=False, centre_points=None, axial_centre_points=None, replicates=1,
              blocks=1, randomize=True, seed=None):
    """Central composite (kind "ccd") or Box-Behnken ("bbd") design.

    centre_points: per block (unblocked: per replicate); for a blocked central
    composite design, per cube block, with axial_centre_points in the axial block.
    None: automatic."""
    k = len(factors)
    if kind not in RSM_KINDS:
        raise DesignError(f"Unknown design {kind!r}.")
    text = [f.name for f in factors if not f.numeric]
    if text:
        raise DesignError("Response surface designs need numeric levels; "
                          + ", ".join(text) + " has text levels.")
    if replicates < 1 or blocks < 1 or (centre_points or 0) < 0 or \
            (axial_centre_points or 0) < 0:
        raise DesignError("Centre points must be 0 or more; replicates and blocks 1 or more.")
    notes, properties = [], []
    if kind == "ccd":
        per_rep = _per_replicate_blocks(blocks, replicates, tuple(CCD_CUBE_BLOCKS),
                                        "Central composite")
        groups, a, cube, centre, axial_centre, block_effects = _ccd_groups(
            k, alpha, custom_alpha, fraction, centre_points, axial_centre_points, per_rep)
    else:
        groups, centre = _bbd_groups(k, centre_points,
                                     _per_replicate_blocks(blocks, replicates,
                                                           (1, 3) if k in BB_BLOCKS else (1,),
                                                           "Box-Behnken"))
        a, cube, axial_centre, block_effects = None, None, None, []

    runs, std = [], 0
    per_rep = len(groups)
    for rep in range(replicates):
        for j, group in enumerate(groups):
            block = rep * per_rep + j + 1 if blocks > 1 else 1
            for ptype, row in group:
                std += 1
                runs.append(Run(std, block, ptype, list(row)))
    total = len(runs)
    if total > MAX_RUNS:
        raise DesignError(f"The design would have {total} runs; the limit is {MAX_RUNS}.")
    n_blocks = blocks if blocks > 1 else 1

    points = [r.coded for r in runs]
    p = model_terms(k)
    properties.append(("Rotatable", "yes" if is_rotatable(points) else "no"))
    orthogonal = blocks_orthogonal(runs)
    if orthogonal is not None:
        properties.append(("Blocks orthogonal to the quadratic model",
                           "yes" if orthogonal else "no"))
    properties.append(("Terms in the full quadratic model", p))
    residual = total - p - (n_blocks - 1)
    if residual < 0:
        raise DesignError(f"{total} runs are too few to fit the {p} terms of the quadratic "
                          "model; add centre points or replicates.")
    properties.append(("Residual degrees of freedom", residual))
    distinct = Counter((r.block, tuple(r.coded)) for r in runs)
    pure = sum(n - 1 for n in distinct.values())
    properties.append(("Pure-error degrees of freedom (repeated points)", pure))
    if pure == 0:
        notes.append("No point is repeated, so lack of fit cannot be tested; add centre "
                     "points.")
    if kind == "ccd" and a > 1 and not inscribed:
        negative = [f.name for f in factors if min(f.low, f.high) >= 0 > f.actual(-a)]
        if negative:
            notes.append("At -alpha the setting of " + ", ".join(negative) + " is negative, "
                         "which is probably impossible. Choose an inscribed design, the "
                         "face-centred alpha, or wider Low-High levels.")
        notes.append("The axial points (coded +-alpha) lie outside the Low-High range; check "
                     "that these settings are possible, or choose an inscribed design.")

    runs, seed = _run_order(runs, n_blocks, randomize, seed)
    base = sum(1 for g in groups for ptype, _ in g if ptype != CENTRE)
    return Design(kind, factors, base, runs, replicates, blocks, centre, cube, block_effects,
                  seed,
                  notes, a, a if inscribed and kind == "ccd" else 1.0, axial_centre,
                  properties)


def _ccd_groups(k, alpha, custom_alpha, fraction, centre_points, axial_centre_points,
                per_rep):
    """Blocks of one replicate as lists of (point type, coded row)."""
    if alpha not in ALPHAS:
        raise DesignError(f"Unknown axial distance {alpha!r}.")
    cube = ccd_cube(k, fraction)
    matrix = cube.matrix()
    n_cube = len(matrix)
    cube_blocks = CCD_CUBE_BLOCKS[per_rep]
    rot = rotatable_alpha(n_cube)
    if alpha == "custom" and (custom_alpha is None or custom_alpha <= 0):
        raise DesignError("Enter a positive axial distance alpha.")
    if alpha == "orthogonal" and per_rep == 1:
        raise DesignError("Orthogonal blocking needs blocks: choose 2 or more blocks per "
                          "replicate, or another axial distance.")

    if per_rep == 1:
        centre = uniform_precision_centre(k, n_cube) if centre_points is None else centre_points
        axial_centre = None
    else:
        # automatic: 2 per block, so every block has its own repeated centre points
        centre = BLOCK_CENTRE if centre_points is None else centre_points
        axial_centre = BLOCK_CENTRE if axial_centre_points is None else axial_centre_points
    a = {"rotatable": rot, "face": 1.0, "custom": custom_alpha}.get(alpha)
    if alpha == "orthogonal":
        a = orthogonal_alpha(k, n_cube, centre * cube_blocks, axial_centre)

    axial = [[(s * a if j == i else 0.0) for j in range(k)] for i in range(k) for s in (-1, 1)]
    if per_rep == 1:
        return ([[(CUBE, r) for r in matrix] + [(AXIAL, r) for r in axial] +
                 [(CENTRE, [0] * k)] * centre], a, cube, centre, None, [])

    gens = block_generators(cube, cube_blocks) if cube_blocks > 1 else []
    if gens and _pattern(cube, subgroup(gens))[1]:
        raise DesignError(f"Splitting the {n_cube} factorial runs into {cube_blocks} blocks "
                          "would confound two-factor interactions with blocks. Use fewer "
                          "blocks.")
    label = {}
    for r in matrix:
        label.setdefault(tuple(math.prod(r[i] for i in range(k) if e >> i & 1) for e in gens),
                         len(label))
    groups = [[] for _ in range(cube_blocks)]
    for r in matrix:
        groups[label[tuple(math.prod(r[i] for i in range(k) if e >> i & 1)
                           for e in gens)]].append((CUBE, r))
    groups = [g + [(CENTRE, [0] * k)] * centre for g in groups]
    groups.append([(AXIAL, r) for r in axial] + [(CENTRE, [0] * k)] * axial_centre)
    return groups, a, cube, centre, axial_centre, gens


def _bbd_groups(k, centre_points, per_rep):
    sets = box_behnken_sets(k)
    if per_rep == 1:
        centre = BB_CENTRE[k] if centre_points is None else centre_points
        rows = [r for s in sets for r in _subfactorial(k, s)]
        return [[(EDGE, r) for r in rows] + [(CENTRE, [0] * k)] * centre], centre
    centre = BLOCK_CENTRE if centre_points is None else centre_points
    groups = []
    for block_sets in BB_BLOCKS[k]:
        rows = [r for s in block_sets for r in _subfactorial(k, s)]
        groups.append([(EDGE, r) for r in rows] + [(CENTRE, [0] * k)] * centre)
    return groups, centre
