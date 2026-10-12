<p align="center">
  <img src="LibreStats_LOGO.svg" alt="LibreStats logo" width="200">
</p>

<h1 align="center">LibreStats</h1>

<p align="center"><strong>A statistics package for LibreOffice Calc, written in pure Python.</strong></p>

LibreStats adds a **LibreStats** menu to Calc with dialog-driven statistical tools, in the spirit of Excel statistics add-ins. Each analysis writes its results to a new sheet: a readable report with the parameters used, the method references and notes on assumptions. You can rerun it, check it and cite it.

The first module is **Data Preparation**: outlier detection, missing-data diagnostics and imputation, transformations and scaling, compositional log-ratios, and coding and recoding of variables.

> **Status:** early development (v0.3.0). Developed and tested on Linux with LibreOffice 26.8. Not yet tested on Windows or macOS.

---

## Features

| Area | Tools |
|---|---|
| **Outliers** | Tukey fences (IQR), robust z-scores with median/MAD (Leys et al., 2013), Grubbs test, Dixon test, Mahalanobis distance |
| **Missing data** | Missing-value summary and patterns, Little's MCAR test; imputation by mean, median, k-nearest neighbours, or multiple imputation (MICE: predictive mean matching or Bayesian linear regression) |
| **Transformations** | log (ln, log10, log2) and square root with offsets for zeros, Box-Cox and Yeo-Johnson (λ by maximum likelihood with 95% CI), logit |
| **Scaling** | Centring, autoscaling (z-score), range 0–1, Pareto (van den Berg et al., 2006) |
| **Compositional data** | CLR and ILR (pivot coordinates) log-ratio transforms, with optional zero replacement (Aitchison, 1982) |
| **Coding** | Dummy variables (one-hot), presence/absence, discretisation into bins (equal width, quantiles, optimal 1-D k-means), recoding and merging categories, ranks, normal scores (Blom, van der Waerden, Tukey, Rankit) |

Across all tools:

- **No dependencies.** Pure Python using only the standard library, so it runs on the Python bundled with LibreOffice (no numpy or scipy needed).
- **Verified results.** The methods are checked against R (`naniar`, `outliers`), scipy, scikit-learn, pandas and scikit-bio. See [Verification](#verification).
- **One undo step per analysis.** Optionally, flagged or imputed cells are highlighted in your data.
- **Scriptable.** Every analysis can also be run without dialogs from a macro (see [documentation.md](documentation.md#automation-macros)).

See **[documentation.md](documentation.md)** for the full user guide: every dialog option, formula and output.

---

## Requirements

- **LibreOffice 26.8 or newer** (enforced by the extension)
- Nothing else. The extension uses LibreOffice's built-in Python.

## Installation

### From a built `.oxt`

1. In LibreOffice, open **Tools → Extension Manager… → Add…** and choose `librestats.oxt`,
   or from a terminal: `unopkg add librestats.oxt`
2. Restart LibreOffice.
3. Open a spreadsheet. The **LibreStats** menu appears in the menu bar (Calc only).

### Building from source

```bash
git clone <repository-url>
cd project_python_calc/librestats
./build.sh                       # creates build/librestats.oxt
unopkg add build/librestats.oxt  # e.g. /opt/libreoffice26.8/program/unopkg
```

To update an installed copy:

```bash
unopkg remove org.librestats && ./build.sh && unopkg add build/librestats.oxt
```

Restart LibreOffice after every update, because it caches Python modules.

## Quick start

1. Put your data in columns, with variable names in the first row.
2. Select the data. If you select a single cell, the range expands to the surrounding data block.
3. Choose, for example, **LibreStats → Data Preparation → Outliers…**
4. Pick the methods and click **OK**. The results appear on a new sheet.

---

## Project layout

```
project_python_calc/
├── README.md
├── documentation.md                # user guide
├── libreoffice-dev-setup.md        # development environment setup (Distrobox, LO 26.8 + SDK)
└── librestats/
    ├── build.sh                    # zips src/ into build/librestats.oxt
    ├── src/                        # contents of the extension
    │   ├── META-INF/manifest.xml
    │   ├── description.xml         # id org.librestats, version, LO >= 26.8
    │   ├── Addons.xcu              # LibreStats menu
    │   ├── librestats_ui.py        # UNO component: menu commands and macro interface
    │   └── pythonpath/librestats/  # the package (on sys.path automatically)
    │       ├── distributions.py    # t, chi-square, normal distributions
    │       ├── linalg.py           # Cholesky and related matrix routines
    │       ├── outliers.py         # outlier methods
    │       ├── dixon_table.py      # generated Dixon critical values
    │       ├── missing.py          # EM, Little's test, imputation
    │       ├── transform.py        # transformations, scaling, ranks, normal scores
    │       ├── compositional.py    # CLR, ILR, zero replacement
    │       ├── coding.py           # dummies, bins, recoding
    │       ├── *_calc.py           # read ranges, run analyses, write report sheets
    │       ├── *_dialog.py, dialogs.py   # dialogs
    │       └── calc.py             # Calc helpers
    ├── tests/                      # unit tests and headless integration test
    └── tools/                      # Dixon table generator, Little's test cross-check
```

The statistical modules (`distributions`, `linalg`, `outliers`, `missing`, `transform`, `compositional`, `coding`) contain no LibreOffice code. You can import and test them with any Python 3.

---

## Development

The development environment is described in [libreoffice-dev-setup.md](libreoffice-dev-setup.md).

### Tests

Unit tests compare the pure-Python methods with reference libraries. Those libraries are needed **only for testing**, in a separate virtual environment:

```bash
python3 -m venv .venv-dev
.venv-dev/bin/pip install numpy scipy scikit-learn pandas scikit-bio
cd librestats
../.venv-dev/bin/python -m unittest discover -s tests -p "test_[mdotc]*.py"
```

The integration test starts a headless LibreOffice, runs every tool through the installed extension and checks the sheets it produces:

```bash
cd librestats
/opt/libreoffice26.8/program/python tests/test_extension.py
```

The test can't display dialogs or check that undo removes inserted sheets (headless LibreOffice limitations). Check these by hand in the normal user interface.

### Regenerating the Dixon table

```bash
cd librestats
../.venv-dev/bin/python tools/make_dixon_table.py 20000000
```

### Verification

| Method | Reference |
|---|---|
| t, χ², normal distributions | scipy, relative error < 10⁻⁹ |
| Grubbs critical values | published tables (e.g. G = 2.290 for n = 10, α = 0.05), and an empirical rejection rate of 5% on outlier-free normal data |
| Dixon critical values | own simulation (20 million samples per n); within 0.005 of published tables |
| Mahalanobis distance | scipy |
| Little's MCAR test | R `naniar::mcar_test(airquality)`: χ² = 35.106, df = 14, p = 0.00142 |
| k-NN imputation | scikit-learn `KNNImputer` (identical values) |
| MICE | simulation: recovers the mean, SD and correlation under MAR missingness |
| Box-Cox, Yeo-Johnson, logit, scaling, ranks | scipy |
| CLR, ILR | scikit-bio |
| Equal-width / quantile bins, dummies | pandas |
| 1-D k-means | exhaustive dynamic programming; never worse than scikit-learn `KMeans` |

---

## Roadmap

- Further modules beyond data preparation (Monte Carlo risk analysis is planned)
- Pooling of multiply-imputed analyses (Rubin's rules)
- Testing on Windows and macOS

## References

The main method references are listed in [documentation.md](documentation.md#references).

## License

Not yet chosen.
