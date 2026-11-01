"""Run the post hoc tests on a Calc range and write the report sheet.

Groups are read as for the ANOVA: one group per column, or one column of
values split by a group range.
"""

from dataclasses import dataclass

from . import calc
from . import posthoc as ph
from .assumptions import _median
from .assumptions_calc import _conclusion, _intro, _notes, _write, read_samples
from .compare_calc import ALTERNATIVE_SIGNS, _anova_samples, _bound, _confidence, _descriptives

REFERENCES = {
    "tukey": ["Tukey, J. W. (1953). The problem of multiple comparisons. Unpublished manuscript, "
              "Princeton University.",
              "Kramer, C. Y. (1956). Biometrics, 12(3), 307-310."],
    "games": ["Games, P. A., & Howell, J. F. (1976). Journal of Educational Statistics, 1(2), "
              "113-125."],
    "dunnett": ["Dunnett, C. W. (1955). Journal of the American Statistical Association, "
                "50(272), 1096-1121.",
                "Dunnett, C. W. (1964). Biometrics, 20(3), 482-491."],
    "dunn": ["Kruskal, W. H., & Wallis, W. A. (1952). Journal of the American Statistical "
             "Association, 47(260), 583-621.",
             "Dunn, O. J. (1964). Technometrics, 6(3), 241-252."],
    "holm": ["Holm, S. (1979). Scandinavian Journal of Statistics, 6(2), 65-70."],
    "sidak": ["Šidák, Z. (1967). Journal of the American Statistical Association, 62(318), "
              "626-633."],
    "bh": ["Benjamini, Y., & Hochberg, Y. (1995). Journal of the Royal Statistical Society B, "
           "57(1), 289-300."],
    "bonferroni": [],
}

ADJUST_NAMES = {"bonferroni": "Bonferroni", "sidak": "Šidák", "holm": "Holm",
                "bh": "Benjamini-Hochberg (false discovery rate)"}


@dataclass
class PosthocOptions:
    range_text: str = ""
    has_labels: bool = True
    group_range: str = ""
    alpha: float = 0.05
    control: str = ""               # Dunnett: control group name; empty = first group
    alternative: str = "two-sided"  # Dunnett: two-sided | greater | less
    adjust: str = "holm"            # Dunn: bonferroni | sidak | holm | bh


def _read(doc, opts):
    samples = read_samples(doc, opts.range_text, opts.has_labels, opts.group_range)
    return samples, _anova_samples(samples)


def _settings(opts, family):
    return [("Comparisons", family),
            ("Family-wise confidence level", _confidence(opts.alpha)),
            ("Significance level", opts.alpha)]


def _params(report, params):
    for label, value in params:
        report.add(label, value)
    report.blank()


def _pair_conclusion(p, alpha):
    return _conclusion(p, alpha, "means differ", "no evidence of a difference")


def run_tukey(doc, opts):
    samples, items = _read(doc, opts)
    res = ph.tukey_hsd([s.values for s in items], opts.alpha)
    k = len(items)
    report = calc.Report()
    _intro(report, "Tukey's HSD (Tukey-Kramer for unequal group sizes)", samples,
           _settings(opts, f"all pairs ({k * (k - 1) // 2})"))
    _descriptives(report, "Group", [(s.name, s.values) for s in items])
    _params(report, res.params)
    conf = _confidence(opts.alpha)
    report.add("Group 1", "Group 2", "Difference", "SE", "q", "p (adjusted)", f"{conf} CI lower",
               f"{conf} CI upper", "Conclusion", bold=True)
    for c in res.comparisons:
        report.add(items[c.first].name, items[c.second].name, c.estimate, c.se, c.statistic, c.p,
                   c.ci[0], c.ci[1], _pair_conclusion(c.p, opts.alpha))
    notes = ["Difference = mean(Group 1) - mean(Group 2). p-values and confidence intervals are "
             "adjusted for all pairwise comparisons (studentized range distribution, "
             "q = |difference| / (SE / √2)).",
             "Tukey's HSD uses the pooled variance and assumes equal variances; with unequal "
             "variances use Games-Howell. With unequal group sizes (Tukey-Kramer) it is "
             "slightly conservative.",
             "Tukey's HSD controls the family-wise error rate by itself; it does not require a "
             "significant ANOVA first."]
    _notes(report, samples, notes, REFERENCES["tukey"])
    return _write(doc, report, "Tukey HSD", "Tukey's HSD")


def run_games_howell(doc, opts):
    samples, items = _read(doc, opts)
    res = ph.games_howell([s.values for s in items], opts.alpha)
    k = len(items)
    report = calc.Report()
    _intro(report, "Games-Howell test (unequal variances)", samples,
           _settings(opts, f"all pairs ({k * (k - 1) // 2})"))
    _descriptives(report, "Group", [(s.name, s.values) for s in items])
    conf = _confidence(opts.alpha)
    report.add("Group 1", "Group 2", "Difference", "SE", "t", "df", "p (adjusted)",
               f"{conf} CI lower", f"{conf} CI upper", "Conclusion", bold=True)
    for c in res.comparisons:
        report.add(items[c.first].name, items[c.second].name, c.estimate, c.se, c.statistic, c.df,
                   c.p, c.ci[0], c.ci[1], _pair_conclusion(c.p, opts.alpha))
    notes = ["Difference = mean(Group 1) - mean(Group 2). Each pair uses its own variances and "
             "Welch degrees of freedom; p-values and intervals use the studentized range "
             "distribution for all pairs.",
             "Games-Howell does not assume equal variances or equal group sizes. It can be "
             "slightly liberal when groups are very small (fewer than about 6 values)."]
    _notes(report, samples, notes, REFERENCES["games"])
    return _write(doc, report, "Games-Howell", "Games-Howell test")


def _control_index(items, control):
    name = control.strip()
    if not name:
        return 0
    for i, s in enumerate(items):
        if s.name == name:
            return i
    names = ", ".join(s.name for s in items)
    raise calc.RangeError(f"The control group '{name}' is not one of the groups ({names}).")


def run_dunnett(doc, opts):
    samples, items = _read(doc, opts)
    control = _control_index(items, opts.control)
    res = ph.dunnett([s.values for s in items], control, opts.alternative, opts.alpha)
    ctrl = items[control].name
    sign = ALTERNATIVE_SIGNS[opts.alternative]
    report = calc.Report()
    _intro(report, "Dunnett's test (each group against a control)", samples,
           [("Control group", ctrl),
            ("Alternative hypothesis", f"mean(group) {sign} mean({ctrl})")]
           + _settings(opts, f"each group against {ctrl} ({len(items) - 1})"))
    _descriptives(report, "Group", [(s.name, s.values) for s in items])
    _params(report, res.params)
    conf = _confidence(opts.alpha)
    report.add("Group", "Difference", "SE", "t", "p (adjusted)", f"{conf} CI lower",
               f"{conf} CI upper", "Conclusion", bold=True)
    for c in res.comparisons:
        report.add(items[c.first].name, c.estimate, c.se, c.statistic, c.p, _bound(c.ci[0]),
                   _bound(c.ci[1]),
                   _conclusion(c.p, opts.alpha, f"mean {sign} control",
                               "no evidence of a difference"))
    notes = [f"Difference = mean(group) - mean({ctrl}). p-values and confidence intervals are "
             "adjusted for the comparisons with the control, using their joint multivariate t "
             "distribution.",
             "Dunnett's test uses the pooled variance of all groups and assumes equal "
             "variances. It is more powerful than comparing all pairs when only comparisons "
             "with the control matter."]
    _notes(report, samples, notes, REFERENCES["dunnett"])
    return _write(doc, report, "Dunnett", "Dunnett's test")


def run_dunn(doc, opts):
    samples, items = _read(doc, opts)
    groups = [s.values for s in items]
    h, df, p = ph.kruskal_wallis(groups)
    res = ph.dunn(groups, opts.adjust)
    k = len(items)
    report = calc.Report()
    _intro(report, "Dunn's test (after Kruskal-Wallis)", samples,
           [("Comparisons", f"all pairs ({k * (k - 1) // 2})"),
            ("p-value adjustment", ADJUST_NAMES[opts.adjust]),
            ("Significance level", opts.alpha)])
    report.add("Group", "n", "Median", "Mean rank", bold=True)
    for s, r in zip(items, res.mean_ranks):
        report.add(s.name, len(s.values), _median(s.values), r)
    report.blank()
    report.add("Kruskal-Wallis test", bold=True)
    report.add("H", "df", "p-value", "Conclusion", bold=True)
    report.add(h, df, p, _conclusion(p, opts.alpha, "distributions differ",
                                     "no evidence of a difference"))
    report.blank()
    report.add("Group 1", "Group 2", "Mean rank difference", "SE", "z", "p (unadjusted)",
               "p (adjusted)", "Conclusion", bold=True)
    for c in res.comparisons:
        report.add(items[c.first].name, items[c.second].name, c.estimate, c.se, c.statistic,
                   c.p_raw, c.p, _conclusion(c.p, opts.alpha, "differ",
                                             "no evidence of a difference"))
    notes = ["Ranks are taken over all groups together (ties get their average rank), and the "
             "standard errors include the tie correction.",
             "Kruskal-Wallis and Dunn's test compare distributions (stochastic dominance); they "
             "compare medians only if the groups have distributions of the same shape.",
             {"bonferroni": "Bonferroni is Dunn's (1964) original adjustment; Holm's is "
                            "uniformly more powerful with the same family-wise error control.",
              "sidak": "Šidák's adjustment is exact for independent tests and slightly more "
                       "powerful than Bonferroni.",
              "holm": "Holm's step-down adjustment controls the family-wise error rate and is "
                      "uniformly more powerful than Bonferroni.",
              "bh": "Benjamini-Hochberg controls the false discovery rate (the expected "
                    "proportion of false positives among the rejections), not the family-wise "
                    "error rate."}[opts.adjust]]
    _notes(report, samples, notes, REFERENCES["dunn"] + REFERENCES[opts.adjust])
    return _write(doc, report, "Dunn", "Dunn's test")


RUNNERS = {"posthocTukey": run_tukey, "posthocGamesHowell": run_games_howell,
           "posthocDunnett": run_dunnett, "posthocDunn": run_dunn}
