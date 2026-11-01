"""Create a design sheet: screening or response surface design.

The factors are read from a table on the sheet (Name, Unit, Low, High). The
"Design" sheet lists the runs in run order with coded and actual levels and
empty response columns, and to the right a summary: the factors, design
generators, defining relation, aliases, the effects confounded with blocks and,
for response surface designs, the axial distance, rotatability, orthogonal
blocking and the degrees of freedom of the quadratic model.
"""

from dataclasses import dataclass

import uno

from . import calc
from . import doe

RESPONSE_COLOR = 0xFFF5CE     # light yellow: cells to fill in
LEVEL_FORMAT = "0.####"       # coded levels such as alpha = 1.68179... shown as 1.6818
REFERENCES = {
    "full": [],
    "fractional": ["Box, G. E. P., & Hunter, J. S. (1961). Technometrics, 3(3), 311-351.",
                   "Fries, A., & Hunter, W. G. (1980). Technometrics, 22(4), 601-608."],
    "pb": ["Plackett, R. L., & Burman, J. P. (1946). Biometrika, 33(4), 305-325."],
    "ccd": ["Box, G. E. P., & Wilson, K. B. (1951). Journal of the Royal Statistical Society "
            "B, 13(1), 1-45.",
            "Box, G. E. P., & Hunter, J. S. (1957). Annals of Mathematical Statistics, 28(1), "
            "195-241."],
    "bbd": ["Box, G. E. P., & Behnken, D. W. (1960). Technometrics, 2(4), 455-475."],
}
POINT_TYPES = {
    "ccd": "Point type: 1 = factorial (cube) point, -1 = axial point, 0 = centre point.",
    "bbd": "Point type: 2 = edge midpoint, 0 = centre point.",
}
ALPHA_NAMES = {"rotatable": "rotatable", "orthogonal": "orthogonal blocking",
               "face": "face-centred", "custom": "custom"}


@dataclass
class DesignOptions:
    range_text: str = ""          # factor table: Name, Unit, Low, High
    has_labels: bool = True
    kind: str = "full"            # full | fractional | pb | ccd | bbd
    runs: int = None              # fractional / Plackett-Burman; None: automatic
    centre_points: int = 0        # per block (central composite in blocks: per cube block);
                                  # response surface designs: None = automatic
    axial_centre_points: int = None   # central composite in blocks; None = automatic
    alpha: str = "rotatable"      # central composite: rotatable | orthogonal | face | custom
    custom_alpha: float = None
    fractional_cube: bool = True  # central composite, 5+ factors: resolution V fraction
    inscribed: bool = False       # central composite: axial points at Low / High
    replicates: int = 1
    blocks: int = 1
    randomize: bool = True
    seed: int = None              # None: a new random seed (reported)
    responses: str = "Response"   # comma-separated names


def response_names(text):
    names = [n.strip() for n in text.split(",") if n.strip()]
    return names or ["Response"]


def design_table(design, responses):
    header = ["Std order", "Run order", "Point type", "Block"]
    header += [f"{f.name} (coded)" for f in design.factors]
    header += [f.label() for f in design.factors] + responses
    rows = [header]
    for r in design.runs:
        rows.append([r.std, r.order, r.ptype, r.block] + list(r.coded) +
                    [f.setting(c, design.scale) for f, c in zip(design.factors, r.coded)] +
                    [""] * len(responses))
    return rows


def summary(design, alpha_kind="rotatable"):
    """(rows whose width sets the column widths, long rows that may overflow, index of
    the factor table header in the first)."""
    k = len(design.factors)
    total = len(design.runs)
    fitted = [["Design", design.title()],
              ["Factors", k],
              ["Runs", total],
              ["Design points (without centre points)", design.base_runs],
              ["Replicates", design.replicates],
              ["Blocks", design.blocks]]
    if design.axial_centre_points is not None:
        fitted += [["Centre points per cube block", design.centre_points],
                   ["Centre points in the axial block", design.axial_centre_points]]
    else:
        fitted.append(["Centre points per block", design.centre_points])
    if design.alpha is not None:
        fitted.append(["Axial distance alpha (coded)",
                       f"{design.alpha:.4f} ({ALPHA_NAMES[alpha_kind]}"
                       + (", inscribed: +-alpha at Low and High)" if design.scale != 1 else ")")])
    fitted += [list(p) for p in design.properties]
    fitted.append(["Run order", f"randomised within blocks, seed {design.seed}"
                   if design.seed is not None else "standard order (not randomised)"])
    fitted.append([])
    header = len(fitted)
    if design.kind == "ccd":
        a = design.alpha
        levels = (-a, -1, 0, 1, a)
        fitted.append(["Factor", "Name", "Unit", "-alpha", "Low (-1)", "Centre (0)",
                       "High (+1)", "+alpha"])
    elif design.kind == "bbd":
        levels = (-1, 0, 1)
        fitted.append(["Factor", "Name", "Unit", "Low (-1)", "Centre (0)", "High (+1)"])
    else:
        levels = None
        fitted.append(["Factor", "Name", "Unit", "Low (-1)", "High (+1)"])
    for i, f in enumerate(design.factors):
        values = [f.low, f.high] if levels is None else \
            [f.setting(c, design.scale) for c in levels]
        fitted.append([doe.LETTERS[i], f.name, f.unit] + values)

    long = [[], [POINT_TYPES.get(design.kind, "Point type: 1 = factorial point, "
                                                "0 = centre point.")]]
    if design.fraction is not None and design.fraction.words:
        if design.kind == "ccd":
            long += [[], [f"Factorial part: 2^({k}-{k - design.fraction.runs.bit_length() + 1})"
                          f", resolution {doe.roman(design.fraction.resolution)}"]]
        long += [[], ["Design generators"]] + [[g] for g in design.generator_text()]
        long += [[], ["Defining relation"], [design.defining_relation()]]
        long += [[], [f"Aliases (main effects and two-factor interactions; interactions up "
                      f"to order {doe.ALIAS_ORDER} shown)"]]
        long += [[c] for c in design.alias_chains()]
    if design.block_effects:
        title = "Confounded with the cube blocks" if design.kind == "ccd" else \
            "Confounded with blocks"
        long += [[], [title]] + [[e] for e in design.block_text()]
    if design.notes:
        long += [[]] + [["Note: " + n] for n in design.notes]
    refs = REFERENCES[design.kind]
    if refs:
        long += [[], ["References"]] + [[r] for r in refs]
    return fitted, long, header


BOLD_TITLES = {"Design generators", "Defining relation", "Confounded with blocks",
               "Confounded with the cube blocks", "References"}


def run_design(doc, opts):
    cell_range = calc.resolve_range(doc, opts.range_text)
    data = calc.used_data(cell_range)
    factors = doe.read_factors(data[1:] if opts.has_labels else data)
    if opts.kind in doe.RSM_KINDS:
        design = doe.build_rsm(factors, opts.kind, opts.alpha, opts.custom_alpha,
                               opts.fractional_cube, opts.inscribed, opts.centre_points,
                               opts.axial_centre_points, opts.replicates, opts.blocks,
                               opts.randomize, opts.seed)
    else:
        design = doe.build_design(factors, opts.kind, opts.runs, opts.centre_points or 0,
                                  opts.replicates, opts.blocks, opts.randomize, opts.seed)
    responses = response_names(opts.responses)
    table = design_table(design, responses)
    fitted, long, factor_header = summary(design, opts.alpha)
    width = len(table[0])
    left = width + 1

    with calc.undo_context(doc, "LibreStats: New design"):
        sheet = calc.new_sheet(doc, "Design")
        data_block = sheet.getCellRangeByPosition(0, 0, width - 1, len(table) - 1)
        data_block.setDataArray(tuple(tuple(calc._cell_value(v) for v in row)
                                      for row in table))
        sheet.getCellRangeByPosition(0, 0, width - 1, 0).CharWeight = calc.BOLD
        k = len(factors)
        formats, locale = doc.getNumberFormats(), uno.createUnoStruct("com.sun.star.lang.Locale")
        key = formats.queryKey(LEVEL_FORMAT, locale, False)
        if key == -1:
            key = formats.addNew(LEVEL_FORMAT, locale)
        sheet.getCellRangeByPosition(4, 1, 3 + k, len(table) - 1).NumberFormat = key
        sheet.getCellRangeByPosition(width - len(responses), 1, width - 1,
                                     len(table) - 1).CellBackColor = RESPONSE_COLOR
        sheet.getCellRangeByPosition(0, 0, width - 1, 0).getColumns().OptimalWidth = True

        # labels and factor table set the widths; the summary values overflow to the right
        _write(sheet, left, 0, [r[:1] for r in fitted[:factor_header]] + fitted[factor_header:])
        sheet.getCellRangeByPosition(left, 0, left, len(fitted) - 1).CharWeight = calc.BOLD
        last = left + len(fitted[factor_header]) - 1
        sheet.getCellRangeByPosition(left, factor_header, last,
                                     factor_header).CharWeight = calc.BOLD
        sheet.getCellRangeByPosition(left, 0, last, 0).getColumns().OptimalWidth = True
        _write(sheet, left + 1, 0, [r[1:] or [""] for r in fitted[:factor_header]])
        # long lines are written after the widths are set, so they overflow to the right
        _write(sheet, left, len(fitted), long)
        for i, row in enumerate(long):
            if row and (row[0] in BOLD_TITLES or row[0].startswith("Aliases")):
                sheet.getCellByPosition(left, len(fitted) + i).CharWeight = calc.BOLD
    controller = doc.getCurrentController()
    controller.setActiveSheet(sheet)
    controller.freezeAtPosition(0, 1)      # keep the header row visible
    return sheet.getName()


def _write(sheet, col, row, rows):
    if not rows:
        return
    w = max(len(r) for r in rows) or 1
    sheet.getCellRangeByPosition(col, row, col + w - 1, row + len(rows) - 1).setDataArray(
        tuple(tuple(calc._cell_value(v) for v in list(r) + [""] * (w - len(r)))
              for r in rows))


# menu commands; each opens the dialog with its design type selected
RUNNERS = {"doeFull": run_design, "doeFractional": run_design,
           "doePlackettBurman": run_design, "doeCCD": run_design, "doeBoxBehnken": run_design}
