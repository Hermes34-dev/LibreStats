# LibreOffice 26.8 Development Setup on Bazzite

*Setup log and reproducible guide for **LibreStats**, a statistics package for LibreOffice Calc (an extension written in Python).*

---

## 1. Why we did this

The goal is to build **LibreStats**, a **statistics package for LibreOffice Calc**, delivered as an extension (an `.oxt` file) written in Python. Monte Carlo risk analysis (in the style of @RISK) is **one module** of the package, not its whole scope. The full list of modules is still to be decided.

The package will add:

- **Custom cell functions** that behave like built-in Calc functions, e.g. `=RISKNORMAL(10;2)` for Monte Carlo inputs, plus statistical functions for the other modules.
- **A UI** (menus, toolbar buttons, dialogs) to run analyses, such as a Monte Carlo simulation, and show the results.

To do that, I need:

| Need | Why |
|---|---|
| **LibreOffice 26.8** (the newest release) | I'm only targeting the latest version, so I develop and test against it. |
| **The LibreOffice SDK** | It contains `unoidl-write`, the tool that compiles the IDL file (the definition of my custom cell functions) into a `.rdb` file. Without it, Calc can't see custom functions. |
| **An isolated, controlled environment** | Bazzite is an immutable system, so I can't install RPMs on it directly. A container keeps LibreOffice, the SDK and their dependencies separate from the main system, and pinned to an exact version. |

**Why a container rather than the Flatpak:** the Flatpak version of LibreOffice exists, but there is **no Flatpak version of the SDK**. Putting both LibreOffice and the SDK in one Distrobox container keeps everything in one place.

**Why the official TDF packages rather than Fedora's own `libreoffice`:** the official packages let me choose the exact version (26.8.x). They install into their own folder (`/opt/libreoffice26.8/`) and never update unless I do it deliberately.

---

## 2. What we've done so far

1. Learned how a LibreOffice extension is structured (see section 4).
2. Created a Fedora container called `lo-dev` with Distrobox.
3. Downloaded the official LibreOffice **26.8.1** RPMs and the matching **SDK** from The Document Foundation.
4. Installed both inside the container.
5. Fixed the missing graphics libraries (the `libXinerama.so.1` error). The Fedora container image is minimal and doesn't include them.
6. Verified the install:
   ```
   LibreOffice 26.8.1.1 bf6ef2662401c7d6dce671ab4deb26ee0af3a7c1
   ```

---

## 3. Step-by-step: how to reproduce this setup from scratch

> Commands marked **[Bazzite]** run in a normal Bazzite terminal.
> Commands marked **[container]** run inside `lo-dev`, where the prompt looks like `[OLI@lo-dev ...]$`.

### Step 1: Create the container **[Bazzite]**

```bash
distrobox create --name lo-dev --image registry.fedoraproject.org/fedora:latest
```

This creates a small Fedora system called `lo-dev`. It **shares your home folder** with Bazzite, so files are visible from both sides.

### Step 2: Enter the container **[Bazzite → container]**

```bash
distrobox enter lo-dev
```

The first time takes a few minutes. When the prompt changes to `...@lo-dev`, you're inside.

### Step 3: Find the newest 26.8 version **[browser]**

Open <https://download.documentfoundation.org/libreoffice/stable/> and note the highest folder starting with `26.8`, e.g. `26.8.1/`.

### Step 4: Make a working folder and set the version **[container]**

```bash
mkdir -p ~/lo-install && cd ~/lo-install
VER=26.8.1          # replace with the version found in step 3
```

### Step 5: Download LibreOffice and the SDK **[container]**

```bash
BASE=https://download.documentfoundation.org/libreoffice/stable/$VER/rpm/x86_64
curl -LO $BASE/LibreOffice_${VER}_Linux_x86-64_rpm.tar.gz
curl -LO $BASE/LibreOffice_${VER}_Linux_x86-64_rpm_sdk.tar.gz
ls -lh
```

Check the sizes: the program is a few hundred MB and the SDK is smaller. A file of only a few KB means the version number has a typo.

*(Optional)* Integrity check: run `sha256sum *.tar.gz` and compare with the SHA-256 shown on libreoffice.org's download page.

### Step 6: Unpack **[container]**

```bash
tar -xzf LibreOffice_${VER}_Linux_x86-64_rpm.tar.gz
tar -xzf LibreOffice_${VER}_Linux_x86-64_rpm_sdk.tar.gz
ls
```

The unpacked folders have a build number in the name (e.g. `LibreOffice_26.8.1.1_...`). That's normal.

### Step 7: Install both **[container]**

```bash
sudo dnf install ./LibreOffice_*_rpm/RPMS/*.rpm
sudo dnf install ./LibreOffice_*_rpm_sdk/RPMS/*.rpm
```

The `./` matters: it installs *these local files*, not Fedora's own LibreOffice.
⚠️ **Never** run `sudo dnf install libreoffice` in this container. That would add a second, different LibreOffice.

### Step 8: Install the missing system libraries **[container]**

The Fedora image has no graphics libraries, and the official RPMs don't pull them in automatically:

```bash
sudo dnf install libXinerama libSM libICE libXrandr cairo cups-libs dbus-glib gtk3 nss mesa-libGL
```

If something is *still* missing later, list all missing libraries in one go:

```bash
ldd /opt/libreoffice26.8/program/*.so 2>/dev/null | grep "not found" | sort -u
dnf provides '*/libNAME.so.X'     # find the package for each one, then dnf install it
```

### Step 9: Verify **[container]**

```bash
/opt/libreoffice26.8/program/soffice --version       # → LibreOffice 26.8.1.x
ls /opt/libreoffice26.8/sdk/bin/unoidl-write         # SDK tool present
rpm -qa | grep -i libreoffice                        # only the official packages
```

✅ If all three look right, the setup is complete.

---

## 4. Background: how a Calc extension works (summary)

- An extension is a **`.oxt` file**, which is just a ZIP with a set structure. LibreOffice unpacks and registers it on install.
- Python code talks to LibreOffice through **PyUNO** (the Python bridge to LibreOffice's UNO component API).
- Typical layout:

```
myrisk.oxt
├── META-INF/manifest.xml   # lists every item and its type (folder name must be UPPERCASE)
├── description.xml         # name, version, dependencies, license
├── XRiskFunctions.rdb      # compiled IDL, generated with unoidl-write (SDK)
├── risk_addin.py           # Python implementation
├── CalcAddIns.xcu          # Function Wizard metadata (names, parameter descriptions)
└── Addons.xcu              # menus / toolbars
```

- **Custom cell functions** need three pieces: an IDL file (compiled to `.rdb`), a Python class implementing it, and `CalcAddIns.xcu`.
- **Long-running analyses** (e.g. the Monte Carlo loop) should be separate commands run from a menu, *not* cell functions. For Monte Carlo, the loop is a separate command (a menu button), *not* a cell function. Each iteration samples the inputs, recalculates the sheet (`calculateAll()`) and reads the output cells.

### Design decisions so far

- **Target only LibreOffice 26.8+.** Enforce it in `description.xml`:
  ```xml
  <dependencies>
    <lo:LibreOffice-minimal-version d:name="LibreOffice 26.8" value="26.8"/>
  </dependencies>
  ```
  (with `xmlns:lo="http://libreoffice.org/extensions/description/2011"` and `xmlns:d="http://openoffice.org/extensions/description/2006"` declared on the root element)
- **Avoid numpy/scipy.** The Windows build of LibreOffice bundles its own Python, which doesn't include them. A pure-Python core (`random`, `statistics`, `math`) keeps the extension portable. Anything heavier (e.g. distribution CDFs, matrix routines for regression) has to be written by us.
- **Development happens inside `lo-dev`.** Claude Code also runs inside the container, so it can build and test directly.
- **Test on real Windows** before release. The HP Victus dual-boots Windows 11.

---

## 5. Useful commands (cheat sheet)

| Task | Command **[container]** |
|---|---|
| Enter the container | `distrobox enter lo-dev` *(from Bazzite)* |
| Open LibreOffice Calc | `/opt/libreoffice26.8/program/soffice --calc` |
| Register an extension folder (development) | `/opt/libreoffice26.8/program/unopkg add ~/dev/myrisk/` |
| Remove an extension | `/opt/libreoffice26.8/program/unopkg remove <extension-id>` |
| List installed extensions | `/opt/libreoffice26.8/program/unopkg list` |
| Compile IDL → RDB | `/opt/libreoffice26.8/sdk/bin/unoidl-write /opt/libreoffice26.8/program/types.rdb X.idl X.rdb` |
| Build LibreStats | `librestats/build.sh` → `librestats/build/librestats.oxt` |
| Reinstall LibreStats | `unopkg remove org.librestats && librestats/build.sh && unopkg add librestats/build/librestats.oxt` |
| Unit tests (needs dev venv) | `cd librestats && ../.venv-dev/bin/python -m unittest discover -s tests -p "test_[mod]*.py"` |
| Integration test (headless) | `cd librestats && /opt/libreoffice26.8/program/python tests/test_extension.py` |

> After changing Python code in an extension, **restart LibreOffice**, because it caches Python modules.

---

## 6. Progress

1. ~~Hello-world cell function and test menu~~ (`=STATSHELLO()`, *Stats → Test Dialog…*): proved the build → install → use cycle. Removed when the project was renamed to LibreStats.
2. **Data preparation → Outliers (v0.1.0):** *LibreStats → Data Preparation → Outliers…* (see section 7).
3. **Data preparation → Missing data (v0.2.0):** diagnostics and imputation (see section 7).
4. **Data preparation → Transform and Scale, Coding (v0.3.0)** (see section 7). The menu is now:

```
LibreStats ▸ Data Preparation ▸ Outliers…
                               ▸ Missing Data ▸ Diagnostics… / Impute Missing Values…
                               ▸ Transform and Scale ▸ Transform… / Scale… / Compositional Log-ratios…
                               ▸ Coding ▸ Dummy Variables… / Presence/Absence… / Discretize into Bins…
                                          / Recode Categories… / Ranks and Normal Scores…
```

Next: decide the remaining modules (Monte Carlo is one of them).

---

## 7. LibreStats layout and modules

```
project_python_calc/
├── .venv-dev/                      # dev-only Python with numpy/scipy, for reference tests (not shipped)
└── librestats/
    ├── build.sh                    # zips src/ into build/librestats.oxt
    ├── src/                        # = contents of the .oxt
    │   ├── META-INF/manifest.xml
    │   ├── description.xml         # id org.librestats, version, LO >= 26.8
    │   ├── Addons.xcu              # LibreStats menu (Calc only)
    │   ├── librestats_ui.py        # UNO component org.librestats.Handler (menu + XJob automation)
    │   └── pythonpath/librestats/  # added to sys.path automatically by LibreOffice
    │       ├── distributions.py    # t and chi-square (pure Python, checked against scipy)
    │       ├── dixon_table.py      # GENERATED Dixon critical values
    │       ├── outliers.py         # the methods (pure Python, no UNO)
    │       ├── outliers_calc.py    # reads the range, writes the report sheet, highlights
    │       ├── outliers_dialog.py  # options dialog
    │       ├── missing.py          # patterns, EM, Little's test, mean/median/kNN/MICE (pure Python)
    │       ├── missing_calc.py     # diagnostics report, imputed-data sheets
    │       ├── missing_dialog.py   # diagnostics and imputation dialogs
    │       ├── transform.py        # log/sqrt/Box-Cox/Yeo-Johnson/logit, scaling, ranks, normal scores
    │       ├── compositional.py    # CLR, ILR (pivot), zero replacement
    │       ├── coding.py           # dummies, presence/absence, bins, recoding rules
    │       ├── prep_calc.py        # runs those tools on a range, writes new columns + report
    │       ├── prep_dialog.py      # their dialogs
    │       ├── linalg.py           # Cholesky etc. for symmetric matrices
    │       ├── dialogs.py          # shared dialog-building helpers
    │       └── calc.py             # Calc/UNO helpers
    ├── tests/                      # unit tests (scipy reference) + headless integration test
    └── tools/
        ├── make_dixon_table.py     # regenerates dixon_table.py by simulation (numpy)
        └── check_little_ml.py      # cross-checks Little's test by direct ML optimisation (scipy)
```

**Methods and choices**

| Method | Details |
|---|---|
| Tukey fences | Q1 − k·IQR, Q3 + k·IQR; k = 1.5 / 2.2 (Hoaglin & Iglewicz) / 3.0. Quartiles = Calc `QUARTILE` (type 7). |
| Robust z (MAD) | (x − median) / (1.4826 · MAD); threshold 2 / **2.5** / 3 (Leys et al., 2013). |
| Grubbs | Two-sided, α = 0.10 / 0.05 / 0.01; optional repetition (with masking warning). p = Bonferroni bound, as in R `outliers`. |
| Dixon | Two-sided, n = 3–30; r10 (3–7), r11 (8–10), r21 (11–13), r22 (14–30). Critical values from our own simulation (20 million samples per n), within 0.005 of the published tables. |
| Mahalanobis | Classical mean and covariance, complete rows only; D² vs χ²(p), α = 0.001 by default. Warns when n is too small for any row to reach the cutoff (max D² = (n−1)²/n). |

The output goes to a new sheet called "Outliers". Highlighting the flagged cells in the input data is optional. The whole run is one undo step.

**Missing data**

| Feature | Details |
|---|---|
| What counts as missing | Empty cells, and text in a numeric column (e.g. "NA"). Columns without any numbers are skipped for analysis but copied to the output. |
| Diagnostics | Missing count and % per variable, number of complete rows, a table of missing-data patterns. |
| Little's MCAR test | ML mean and covariance by EM; d² = Σ nⱼ (ȳⱼ − μ̂ⱼ)ᵀ Σ̂ⱼ⁻¹ (ȳⱼ − μ̂ⱼ), df = Σ pⱼ − p (Little, 1988). Matches R `naniar::mcar_test(airquality)` (χ² = 35.106, df = 14, p = 0.00142). |
| Mean / median | Per variable, with a warning that it shrinks variances and correlations. |
| k-nearest neighbours | k = 3 / 5 / 10; Euclidean distance over variables observed in both rows (scikit-learn's `nan_euclidean`), on z-scores by default; donor mean. Matches scikit-learn `KNNImputer`. |
| MICE | Chained equations with Bayesian linear regression (Rubin, 1987, as `mice.impute.norm`) or predictive mean matching (default, 5 donors, as R `mice`); m = 5–50 datasets, 5–20 iterations, seed for reproducibility. Output stacked in long format with "Imputation" and "Row" columns. |

Single imputation can write to a new sheet ("Imputed") or fill the input range in place. MICE always writes to a new sheet ("MICE"). Imputed cells are highlighted in light blue.

**Transform and Scale, Coding**

Each of these tools adds new columns. They go either on a new sheet (a copy of the data plus the new columns, with a report beside it) or directly to the right of the input range; that second option refuses to overwrite non-empty cells. The report lists the parameters used (λ, mean, SD, bin limits, …) so the same transformation can be applied to new data.

| Tool | Details |
|---|---|
| Log | ln / log10 / log2, offset none, +1, + half the smallest positive value, or a custom value. |
| Square root | Same offsets. |
| Box-Cox | λ by maximum likelihood with a 95 % likelihood-ratio CI (or fixed λ); positive values only, optional offset. Matches `scipy.stats.boxcox`. |
| Yeo-Johnson | Any real values; λ by ML with CI. Matches `scipy.stats.yeojohnson`. |
| Logit | Proportions or percentages; exact 0/1 handled by squeezing (Smithson & Verkuilen, 2006). |
| Scaling | Centring, autoscaling (z), range 0–1 (min–max), Pareto (van den Berg et al., 2006). |
| Log-ratios | CLR, or ILR pivot coordinates (Filzmoser et al., 2018); zeros: refuse, or multiplicative replacement with 0.65 × the smallest positive value of the part (Martín-Fernández et al., 2003). Matches scikit-bio. |
| Dummy variables | One 0/1 column per category (numbers or text), optionally leaving out the first category as the reference. |
| Presence/absence | 1 if the value is above a threshold; empty cells as missing or absent. |
| Bins | Equal width or quantiles (right-closed, same as pandas `cut`/`qcut`), or optimal 1-D k-means (Wang & Song, 2011). |
| Recode | Rules such as `1; 2 -> low`, `* -> other`. A button lists the categories in the range as a starting point. |
| Ranks / normal scores | Ranks with average/min/max/dense ties; normal scores by Blom, van der Waerden, Tukey or Rankit. |
Analyses can also be run without the dialog through `org.librestats.Handler`'s `XJob.execute` (see `tests/test_extension.py`).

---

## 8. Reference documentation

- LibreOffice API reference: <https://api.libreoffice.org>
- Extension file format (DevGuide): <https://wiki.openoffice.org/wiki/Documentation/DevGuide/Extensions/File_Format>
- `description.xml` (DevGuide): <https://wiki.openoffice.org/wiki/Documentation/DevGuide/Extensions/description.xml>
- Calc Add-In in Python (worked example): <https://wiki.openoffice.org/wiki/Calc/Add-In/Python_How-To>
- Same tutorial, alternative version: <https://www.biochemfusion.com/doc/Calc_addin_howto.html>
- Modern real-world example with a custom Calc function: <https://github.com/KeithCu/localwriter>
- Extension project generator: <https://github.com/bastien34/mastering_cookiecutter_OOo_extension>
- LibreOffice 26.8 release notes: <https://wiki.documentfoundation.org/ReleaseNotes/26.8>
