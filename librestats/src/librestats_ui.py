"""LibreStats UNO component: handles the LibreStats menu.

Menu items use URLs of the form service:org.librestats.Handler?<command>,
which call trigger(<command>). The same service implements XJob so analyses
can also be run without a dialog (from macros or tests), with named
arguments Command, Document, Range, HasLabels and command-specific options:

  Command="outliers"          Methods="tukey,mad,grubbs,dixon,mahalanobis", TukeyK,
                              MadThreshold, Alpha, GrubbsIterate, MahalanobisAlpha, Highlight
  Command="missing_diagnose"  LittleTest
  Command="impute"            Method="mean|median|knn|mice", K, Standardize,
                              MiceMethod="pmm|norm", M, Iterations, Seed, InPlace, Highlight
  Command="transform", "scale", "logratio", "dummies", "presence", "bins", "recode",
          "ranks"             Beside, plus the lower-case settings listed in prep_calc
  Command="histogram"         Bins="auto|count|width", K, Width, Min, Max, Percent, Title,
                              XTitle, YTitle, Color, Beside
  Command="qq"                GroupRange, Positions="blom|vdw|tukey|rankit",
                              Line="quartiles|meansd|none", Title, XTitle, YTitle, Color
  Command="normality"         GroupRange, Alpha
  Command="variances"         GroupRange, Center="median|mean", Alpha

execute() returns the name of the sheet that was created.
"""

import traceback

import unohelper
from com.sun.star.task import XJob, XJobExecutor

from librestats import (assumptions_calc, assumptions_dialog, calc, missing, missing_calc,
                        missing_dialog, outliers_calc, outliers_dialog, prep_calc, prep_dialog,
                        viz_calc, viz_dialog)

IMPLEMENTATION_NAME = "org.librestats.HandlerImpl"
SERVICE_NAMES = ("org.librestats.Handler",)

# menu command -> (dialog, analysis)
COMMANDS = {
    "outliers": (outliers_dialog.show, outliers_calc.run),
    "missingDiagnose": (missing_dialog.show_diagnose, missing_calc.diagnose),
    "missingImpute": (missing_dialog.show_impute, missing_calc.impute),
}
# transformations and coding: same name for menu command and XJob Command
COMMANDS.update({name: (prep_dialog.DIALOGS[name], prep_calc.RUNNERS[name])
                 for name in prep_calc.RUNNERS})
COMMANDS.update({name: (viz_dialog.DIALOGS[name], viz_calc.RUNNERS[name])
                 for name in viz_calc.RUNNERS})
COMMANDS.update({name: (assumptions_dialog.DIALOGS[name], assumptions_calc.RUNNERS[name])
                 for name in assumptions_calc.RUNNERS})


class LibreStatsHandler(unohelper.Base, XJobExecutor, XJob):
    def __init__(self, ctx):
        self.ctx = ctx

    def _doc_and_window(self):
        desktop = self.ctx.ServiceManager.createInstanceWithContext(
            "com.sun.star.frame.Desktop", self.ctx)
        doc = desktop.getCurrentComponent()
        frame = desktop.getCurrentFrame()
        return doc, frame.getContainerWindow() if frame else None

    # --- menu ---------------------------------------------------------------

    def trigger(self, command):
        doc, window = self._doc_and_window()
        if command not in COMMANDS:
            return
        show, analyse = COMMANDS[command]
        try:
            opts = show(self.ctx, doc, window)
            if opts is not None:
                analyse(doc, opts)
        except (calc.RangeError, missing.MissingDataError, ValueError) as exc:  # incl. method errors
            calc.message_box(self.ctx, window, "LibreStats", str(exc))
        except Exception:
            calc.message_box(self.ctx, window, "LibreStats - unexpected error",
                             traceback.format_exc())

    # --- automation (no dialog) -----------------------------------------------

    def execute(self, arguments):
        args = {a.Name: a.Value for a in arguments}
        command = args.get("Command")
        doc = args.get("Document") or self._doc_and_window()[0]
        common = dict(range_text=args["Range"], has_labels=bool(args.get("HasLabels", True)))

        if command == "outliers":
            methods = {m.strip().lower() for m in str(args.get("Methods", "tukey,mad")).split(",")}
            opts = outliers_calc.OutlierOptions(
                **common,
                tukey="tukey" in methods, tukey_k=float(args.get("TukeyK", 1.5)),
                mad="mad" in methods, mad_threshold=float(args.get("MadThreshold", 2.5)),
                grubbs="grubbs" in methods,
                grubbs_iterate=bool(args.get("GrubbsIterate", False)),
                dixon="dixon" in methods, alpha=float(args.get("Alpha", 0.05)),
                mahalanobis="mahalanobis" in methods,
                mahalanobis_alpha=float(args.get("MahalanobisAlpha", 0.001)),
                highlight=bool(args.get("Highlight", False)))
            return outliers_calc.run(doc, opts)

        if command == "missing_diagnose":
            opts = missing_calc.DiagnoseOptions(**common, little=bool(args.get("LittleTest", True)))
            return missing_calc.diagnose(doc, opts)

        if command == "impute":
            seed = args.get("Seed", 12345)
            opts = missing_calc.ImputeOptions(
                **common,
                method=str(args.get("Method", "mean")),
                k=int(args.get("K", 5)),
                standardize=bool(args.get("Standardize", True)),
                mice_method=str(args.get("MiceMethod", "pmm")),
                m=int(args.get("M", 5)),
                iterations=int(args.get("Iterations", 10)),
                seed=None if seed in (None, "") else int(seed),
                in_place=bool(args.get("InPlace", False)),
                highlight=bool(args.get("Highlight", True)))
            return missing_calc.impute(doc, opts)

        if command in prep_calc.RUNNERS:
            reserved = {"Command", "Document", "Range", "HasLabels", "Beside"}
            opts = prep_calc.PrepOptions(**common, beside=bool(args.get("Beside", False)),
                                         settings={k: v for k, v in args.items()
                                                   if k not in reserved})
            return prep_calc.RUNNERS[command](doc, opts)

        if command == "histogram":
            def number(key):
                v = args.get(key)
                return None if v in (None, "") else float(v)

            opts = viz_calc.HistogramOptions(
                **common, bins=str(args.get("Bins", "auto")), k=int(args.get("K", 10)),
                width=number("Width"), min=number("Min"), max=number("Max"),
                percent=bool(args.get("Percent", False)), title=str(args.get("Title", "")),
                x_title=str(args.get("XTitle", "")), y_title=str(args.get("YTitle", "")),
                color=int(args.get("Color", viz_calc.DEFAULT_COLOR)),
                beside=bool(args.get("Beside", False)))
            return viz_calc.run_histogram(doc, opts)

        if command == "qq":
            opts = viz_calc.QQOptions(
                **common, group_range=str(args.get("GroupRange", "")),
                positions=str(args.get("Positions", "blom")),
                line=str(args.get("Line", "quartiles")), title=str(args.get("Title", "")),
                x_title=str(args.get("XTitle", "")), y_title=str(args.get("YTitle", "")),
                color=int(args.get("Color", viz_calc.DEFAULT_COLOR)))
            return viz_calc.run_qq(doc, opts)

        if command in ("normality", "variances"):
            group = str(args.get("GroupRange", ""))
            alpha = float(args.get("Alpha", 0.05))
            if command == "normality":
                opts = assumptions_calc.NormalityOptions(**common, group_range=group, alpha=alpha)
            else:
                opts = assumptions_calc.VarianceOptions(
                    **common, group_range=group, center=str(args.get("Center", "median")),
                    alpha=alpha)
            return assumptions_calc.RUNNERS[command](doc, opts)

        raise ValueError(f"Unknown LibreStats command: {command!r}")


g_ImplementationHelper = unohelper.ImplementationHelper()
g_ImplementationHelper.addImplementation(LibreStatsHandler, IMPLEMENTATION_NAME, SERVICE_NAMES)
