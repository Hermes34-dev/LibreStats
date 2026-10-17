# LibreStats Documentation

User guide for **LibreStats 0.4.0**, a statistics extension for LibreOffice Calc (26.8 or newer).

## Contents

1. [General concepts](#1-general-concepts)
2. [Outliers](#2-outliers)
3. [Missing data](#3-missing-data)
   - [Diagnostics](#31-diagnostics)
   - [Impute missing values](#32-impute-missing-values)
4. [Transform and scale](#4-transform-and-scale)
   - [Transform](#41-transform)
   - [Scale](#42-scale)
   - [Compositional log-ratios](#43-compositional-log-ratios)
5. [Coding](#5-coding)
   - [Dummy variables](#51-dummy-variables)
   - [Presence/absence](#52-presenceabsence)
   - [Discretize into bins](#53-discretize-into-bins)
   - [Recode categories](#54-recode-categories)
   - [Ranks and normal scores](#55-ranks-and-normal-scores)
6. [Data visualization](#6-data-visualization)
   - [Histogram](#61-histogram)
7. [Automation (macros)](#automation-macros)
8. [Limitations](#8-limitations)
9. [References](#references)

---

## 1. General concepts

### The menu

LibreStats adds a top-level **LibreStats** menu to Calc (it does not appear in Writer, Impress, etc.):

```
LibreStats ▸ Data Preparation ▸ Outliers…
                               ▸ Missing Data ▸ Diagnostics…
                                              ▸ Impute Missing Values…
                               ▸ Transform and Scale ▸ Transform…
                                                     ▸ Scale…
                                                     ▸ Compositional Log-ratios…
                               ▸ Coding ▸ Dummy Variables…
                                        ▸ Presence/Absence…
                                        ▸ Discretize into Bins…
                                        ▸ Recode Categories…
                                        ▸ Ranks and Normal Scores…
           ▸ Data Visualization ▸ Histogram…
```

### Data layout

- **Variables in columns, observations in rows.**
- **Variable names in the first row** (optional). Every dialog has a *First row contains variable names* box, which is ticked automatically when the first row contains text. Without names, columns are called by their letter (e.g. "Column B" or "B").
- **Input range.** Every dialog starts with the current selection. If a single cell is selected, the range expands to the surrounding block of data (like Calc's own *Data → Statistics* tools). You can edit the range in the dialog. A range without a sheet name (e.g. `A1:C20`) refers to the active sheet.
- **Whole columns** (e.g. `A:C`) can be selected; trailing empty rows are ignored.

### Empty and text cells

| Tool | Empty cell | Text in a numeric column |
|---|---|---|
| Outliers | ignored | ignored |
| Missing data | missing | counted as missing (e.g. `NA`), reported |
| Transform, Scale, Log-ratios, Bins, Ranks | output left empty | treated as missing, reported |
| Dummy variables, Recode | output left empty | a category like any other |
| Presence/absence | missing, or absent if chosen | output left empty, reported |
| Histogram | ignored | ignored, reported |

### Output

- Results go to a **new sheet** with a descriptive name ("Outliers", "Missing data", "Transformed", …). If a sheet with that name exists, a number is added ("Outliers 2").
- Reports contain the input range, the settings, the estimated parameters, notes on assumptions, and method references.
- Tools that create new variables (sections 4 and 5) offer two output options:
  - **New sheet:** a copy of the input data with the new columns added on the right, and the report beside it.
  - **Directly to the right of the input range:** the new columns are written next to your data. This is **refused if those cells are not empty**, so nothing is overwritten. The report goes on its own sheet (e.g. "Transformed report").
- Numbers are written at full precision; use Calc's number formats to round the display.

### Undo

Each analysis is **one undo step** (*Edit → Undo "LibreStats: …"*), including the new sheet and any highlighting.

### Highlight colours

| Colour | Meaning |
|---|---|
| Light orange | flagged outlier (optional) |
| Light blue | imputed value (optional, on by default) |

### Errors

If a method cannot be applied to one variable (too few values, zero variance, negative values for a log, …), that variable is reported as *not computed* with the reason, and the other variables are still processed. Problems with the input itself (invalid range, no numbers) are shown in a message box and the dialog reopens.

---

## 2. Outliers

**LibreStats → Data Preparation → Outliers…**

Flags unusual values with up to five methods at once. Univariate methods test **each column separately**; the Mahalanobis distance uses **all columns together**.

### Dialog

| Option | Choices | Default |
|---|---|---|
| Tukey fences (IQR) | k = 1.5 (Tukey) / 2.2 (Hoaglin & Iglewicz) / 3.0 (far out) | on, 1.5 |
| Robust z-scores (median / MAD) | \|z\| > 2.0 (lenient) / 2.5 (recommended) / 3.0 (strict) | on, 2.5 |
| Grubbs test | optional *Repeat until none found* | off |
| Dixon test | 3 to 30 values | off |
| Significance level (Grubbs, Dixon) | 0.10 / 0.05 / 0.01 | 0.05 |
| Mahalanobis distance | α = 0.001 (usual) / 0.01 / 0.05 | off, 0.001 |
| Highlight outliers in the input data | | off |

### Methods

**Tukey fences** (Tukey, 1977). A value is flagged if it lies outside
[Q1 − k·IQR, Q3 + k·IQR], with IQR = Q3 − Q1. Quartiles are computed like Calc's `QUARTILE` (linear interpolation, Hyndman & Fan type 7), so you can reproduce them in the sheet. The score is the distance beyond the nearer quartile in IQR units. Needs at least 3 values. If IQR = 0, every value different from the quartiles is flagged (noted in the report).

**Robust z-scores** (Leys et al., 2013). z = (x − median) / MAD, where MAD = 1.4826 × median(|x − median|). The constant 1.4826 makes MAD estimate the standard deviation under normality. Unlike the mean and SD, the median and MAD are not distorted by the outliers themselves. Thresholds follow Leys et al.: 3 (very conservative), 2.5 (moderately conservative, recommended), 2 (poorly conservative). Cannot be used when at least half the values are identical (MAD = 0).

**Grubbs test** (Grubbs, 1969). Two-sided test of whether the most extreme value is an outlier, assuming the remaining data are normal.
G = max|x − x̄| / s, compared with

  G_crit = ((n − 1)/√n) · √( t² / (n − 2 + t²) ),  t = t(α/(2n), n − 2).

The p-value is the Bonferroni bound 2n·P(T > t) (as in R's `outliers` package). It is exact for small p and shown as 1 when the bound exceeds 1. With *Repeat until none found*, the flagged value is removed and the test repeated; the report shows every step. Repeated testing can miss outliers that mask each other (noted in the report). Needs n ≥ 3.

**Dixon test** (Dixon, 1951). Two-sided test of the most extreme value using gap/range ratios, for 3–30 values. The ratio depends on n, as Dixon recommended:

| n | Ratio (largest value; the smallest is the mirror image) |
|---|---|
| 3–7 | r10 = (x₍ₙ₎ − x₍ₙ₋₁₎) / (x₍ₙ₎ − x₍₁₎) |
| 8–10 | r11 = (x₍ₙ₎ − x₍ₙ₋₁₎) / (x₍ₙ₎ − x₍₂₎) |
| 11–13 | r21 = (x₍ₙ₎ − x₍ₙ₋₂₎) / (x₍ₙ₎ − x₍₂₎) |
| 14–30 | r22 = (x₍ₙ₎ − x₍ₙ₋₂₎) / (x₍ₙ₎ − x₍₃₎) |

Critical values were computed by LibreStats with 20 million simulated normal samples per n (precision about ±0.001). They agree with published tables within 0.005. The p-value is reported as a bracket (e.g. "p < 0.01").

**Mahalanobis distance** (Mahalanobis, 1936). For each row, D² = (x − x̄)ᵀ S⁻¹ (x − x̄), using the sample mean and covariance of all complete rows. Rows with D² above the χ² critical value (df = number of variables) are flagged; α = 0.001 is the usual convention. Needs at least 2 variables and p + 2 complete rows; rows with any missing value are excluded. The report warns when the sample is too small for any row to reach the critical value: within a sample, D² cannot exceed (n − 1)²/n.

### Output: sheet "Outliers"

For each method: a summary table per variable (n, quartiles/fences, median/MAD/limits, …), the list of flagged values with their **cell address**, value, side (low/high) and score, plus notes. Grubbs and Dixon show test statistics, critical values and p-values. Mahalanobis lists the flagged rows with D², p and their values.

---

## 3. Missing data

### 3.1 Diagnostics

**LibreStats → Data Preparation → Missing Data → Diagnostics…**

Option: **Little's MCAR test** (on by default).

#### Output: sheet "Missing data"

- **Missing values by variable:** n observed, number and percentage missing, totals, and the number of complete rows.
- **Missing-data patterns:** each distinct combination of observed (1) and missing (0) variables, with its count and percentage, most frequent first.
- **Little's MCAR test** (Little, 1988): χ², df, p, number of patterns and rows used, with an interpretation.

#### Little's MCAR test

Tests the null hypothesis that data are **missing completely at random** (MCAR). The mean μ and covariance Σ are estimated by maximum likelihood with the EM algorithm (Dempster et al., 1977), using all rows. For each missing-data pattern j with nⱼ rows and observed variables oⱼ:

  d² = Σⱼ nⱼ (ȳⱼ − μ̂ₒⱼ)ᵀ Σ̂ₒⱼ⁻¹ (ȳⱼ − μ̂ₒⱼ),  df = Σⱼ |oⱼ| − p.

Reading the result:

- **p < 0.05:** reject MCAR. Missingness is related to the data (MAR or MNAR), and complete-case analysis or mean imputation can bias results.
- **p ≥ 0.05:** no evidence against MCAR. This **does not prove** MCAR; the test may lack power.
- The test assumes multivariate normality and **cannot distinguish MAR from MNAR**; no test based on the observed data can.
- Rows with every value missing are excluded (reported).

Reproduces R `naniar::mcar_test(airquality)`: χ² = 35.106, df = 14, p = 0.00142.

### 3.2 Impute missing values

**LibreStats → Data Preparation → Missing Data → Impute Missing Values…**

Only columns that contain numbers are imputed. Other columns (IDs, groups) are copied unchanged.

#### Dialog

| Option | Choices | Default |
|---|---|---|
| Method | Mean / Median / k-nearest neighbours / MICE | Mean |
| k (k-NN) | 3 / 5 / 10 | 5 |
| Standardise variables before computing distances (k-NN) | | on |
| MICE method | Predictive mean matching / Bayesian linear regression | PMM |
| Imputations (MICE) | 5 / 10 / 20 / 50 | 5 |
| Iterations (MICE) | 5 / 10 / 20 | 10 |
| Random seed (MICE) | whole number; empty = different each run | 12345 |
| Output | Completed data on a new sheet / Fill the missing cells in the input range | new sheet |
| Highlight imputed cells | | on |

#### Methods

**Mean / median.** Each missing value is replaced by its variable's mean or median. This is simple, but it **shrinks variances and correlations**; use it only when few values are missing (the report says so).

**k-nearest neighbours** (Troyanskaya et al., 2001). For each missing value, the k rows most similar to the incomplete row, among those that have that variable observed, are found. The value is replaced by the mean of their values. The distance is Euclidean over the variables observed in both rows, scaled up for the variables not compared (√(p/q · Σ(dᵢ²)), as scikit-learn's `nan_euclidean`). With *Standardise*, distances use z-scores so that variables measured in large units do not dominate. If no row is comparable, the variable mean is used (reported). Results are identical to scikit-learn's `KNNImputer`.

**MICE: multiple imputation by chained equations** (van Buuren & Groothuis-Oudshoorn, 2011).

1. Missing values start as random draws from the observed values of their variable.
2. Each incomplete variable in turn is regressed on all other variables (using the rows where it is observed). The regression coefficients and residual SD are **drawn from their posterior distribution** (Rubin, 1987), so that the uncertainty of the model is reflected.
3. Missing values are replaced by:
   - **Predictive mean matching** (default): the observed value of one of the 5 cases whose predicted value is closest, chosen at random. Imputed values are always real, observed values (no impossible values such as negative concentrations).
   - **Bayesian linear regression:** the prediction plus random normal noise.
4. Steps 2–3 are repeated for the chosen number of iterations; m independent runs give m completed datasets.

Each variable with missing values needs at least p + 2 observed values (p = number of numeric variables). The same seed always gives the same result.

> **How to use MICE results:** analyse each of the m datasets separately and combine the estimates with Rubin's rules. Do **not** average the imputations into one dataset; that would remove the between-imputation uncertainty that multiple imputation exists to capture.

#### Output

| Method / output option | Result |
|---|---|
| Mean, median, k-NN to a new sheet | sheet **"Imputed"**: completed data at A1 (imputed cells light blue), summary to the right |
| Mean, median, k-NN in place | missing cells of the input range filled (and highlighted); summary on sheet **"Imputation"** |
| MICE (always a new sheet) | sheet **"MICE"**: the m datasets stacked (long format), with columns *Imputation* (1…m) and *Row* (original sheet row); information beside it |

---

## 4. Transform and scale

All tools in this section create new columns (see [Output](#output)) and report the parameters used, so that the same transformation can be applied to new data.

### 4.1 Transform

**LibreStats → Data Preparation → Transform and Scale → Transform…**

| Transformation | New column | Requirement |
|---|---|---|
| Logarithm, base e / 10 / 2 | `ln(X)`, `log10(X)`, `log2(X)`, e.g. `ln(X+1)` with offset | x + offset > 0 |
| Square root | `sqrt(X)` | x + offset ≥ 0 |
| Box-Cox | `BoxCox(X)` | x + offset > 0 |
| Yeo-Johnson | `YeoJohnson(X)` | any real values |
| Logit | `logit(X)` | proportions 0–1 (or percentages 0–100) |

#### Offsets (log, square root, Box-Cox)

Zeros cannot be log-transformed. Choose an offset added to every value first:

- **none**
- **+1**, i.e. log(x + 1)
- **+ half the smallest positive value** of that column (common for concentrations with zeros below the detection limit; computed per column and reported)
- **custom:** any number

#### Box-Cox (Box & Cox, 1964)

  y = (x^λ − 1)/λ for λ ≠ 0,  y = ln x for λ = 0.

Leave *Lambda* empty to estimate λ by **maximum likelihood**: the value that makes the transformed data most nearly normal. The report gives λ and its **95% likelihood-ratio confidence interval**. If the interval contains a simple value (0 = log, 0.5 = square root, 1 = no transformation), you may prefer that value for interpretability and enter it as a fixed λ. λ is searched between −10 and 10; reaching a limit is reported. Results match `scipy.stats.boxcox`.

#### Yeo-Johnson (Yeo & Johnson, 2000)

An extension of Box-Cox that also handles **zero and negative values** without an offset:

| | λ ≠ 0 (or ≠ 2) | λ = 0 (or 2) |
|---|---|---|
| x ≥ 0 | ((x + 1)^λ − 1)/λ | ln(x + 1) |
| x < 0 | −((1 − x)^(2−λ) − 1)/(2 − λ) | −ln(1 − x) |

λ is estimated by maximum likelihood with a 95% CI, or fixed. Results match `scipy.stats.yeojohnson`.

#### Logit

  logit(p) = ln(p / (1 − p)).

- Tick *values are percentages* for data on a 0–100 scale (divided by 100 first).
- The logit of exactly 0 or 1 is infinite. With *squeeze* (default), if any 0 or 1 is present, **all** values of that column are compressed with p′ = (p(n − 1) + 0.5)/n (Smithson & Verkuilen, 2006); this is reported. Without it, such columns are refused.

### 4.2 Scale

**LibreStats → Data Preparation → Transform and Scale → Scale…**

Following van den Berg et al. (2006). The SD uses n − 1.

| Method | Formula | New column | Use |
|---|---|---|---|
| Centring only | x − mean | `centred(X)` | removes offsets, keeps units and variance |
| Autoscaling (z-score) | (x − mean) / SD | `z(X)` | all variables equally important (unit variance) |
| Range (0–1) | (x − min) / (max − min) | `range01(X)` | bounded 0–1 scale |
| Pareto | (x − mean) / √SD | `pareto(X)` | compromise: reduces the weight of large values less than autoscaling, keeps closer to the original structure |

> Note: van den Berg et al. (2006) define *range scaling* as (x − mean)/(max − min). LibreStats' "Range (0–1)" is min–max scaling to the interval 0–1.

The report lists the mean, SD, minimum and maximum used for each variable.

### 4.3 Compositional log-ratios

**LibreStats → Data Preparation → Transform and Scale → Compositional Log-ratios…**

For **compositional data**: parts of a whole, such as fatty-acid profiles in % of total fatty acids, mineral or nutrient composition, or relative abundances. Because the parts sum to a constant, they are not free to vary independently. Standard statistics on raw percentages produce **spurious correlations** and are not scale-invariant (Aitchison, 1982). Log-ratio transforms move the data into ordinary real space, where standard methods are valid.

**Data layout:** each row is one composition (sample), and **every selected column is a part**. Select only the part columns (no ID or group columns). Rows need not sum to exactly 100%; log-ratios do not depend on the total. Rows with a missing part are skipped.

| Transformation | Formula | New columns |
|---|---|---|
| **CLR**, centred log-ratio | clrᵢ = ln(xᵢ / g(x)), g = geometric mean of all D parts | D columns `clr(part)` |
| **ILR**, isometric log-ratio, pivot coordinates | zᵢ = √((D−i)/(D−i+1)) · ln(xᵢ / g(xᵢ₊₁,…,x_D)), i = 1…D−1 | D − 1 columns `ilr1: A vs (B, C, …)`, … |

- **CLR** keeps one value per part (easy to read), but the values in each row sum to 0, so their **covariance matrix is singular**. Avoid CLR for methods that invert it (regression with all parts, MANOVA, Mahalanobis distance).
- **ILR** gives D − 1 coordinates in an orthonormal basis: a full-rank, isometric representation (distances equal Aitchison distances). LibreStats uses **pivot coordinates** (Filzmoser et al., 2018): coordinate 1 contains *all* relative information about the first part versus the others. **Put the part of interest in the first column.** Other software (e.g. R `compositions::ilr`) uses a different basis by default, which gives different numbers but the same geometry.

**Zeros.** Logarithms need positive values. By default, zeros are refused. With *Replace zeros*, each zero in part j is replaced by δⱼ = 0.65 × the smallest positive value of part j (a stand-in for the detection limit). The non-zero parts of that row are shrunk so the row total and their ratios are preserved: multiplicative replacement (Martín-Fernández et al., 2003).

Results match scikit-bio (`clr`, `ilr` with the same basis).

---

## 5. Coding

### 5.1 Dummy variables

**LibreStats → Data Preparation → Coding → Dummy Variables…**

One-hot coding: for each selected column, one **0/1 column per category** (`Group=a`, `Group=b`, …). Categories can be text or numbers and are sorted (numbers first, by value; then text, alphabetically). Empty cells give empty dummies.

*Leave out the first category* produces k − 1 columns. The first category becomes the reference (named in the report), as needed for regression with an intercept, to avoid perfect collinearity.

### 5.2 Presence/absence

**LibreStats → Data Preparation → Coding → Presence/Absence…**

Converts abundances or concentrations to **1 = present** (value above the threshold, default 0) and **0 = absent**. Empty cells stay missing, unless *Empty cells mean absent* is ticked. Text cells are left empty and counted in the report. New column: `X present`.

### 5.3 Discretize into bins

**LibreStats → Data Preparation → Coding → Discretize into Bins…**

Converts a numeric variable into bin numbers 1…k (k = 2–10, default 4). New column: `X bin`.

| Method | Boundaries |
|---|---|
| **Equal width** | k intervals of equal width between the minimum and maximum |
| **Quantiles** | boundaries at the 1/k, 2/k, … quantiles (as Calc `PERCENTILE`), giving about equal counts per bin. If ties make boundaries equal, fewer bins are created (reported). |
| **k-means** | the grouping into k bins with the smallest within-bin sum of squares: natural groups of similar values. Solved **exactly** by dynamic programming (Wang & Song, 2011), not by the random-start algorithm, so the result is optimal and reproducible. |

Equal-width and quantile bins are **right-closed**: (from, to], with the first bin also including its lower limit. This matches Calc's `FREQUENCY` and pandas `cut`/`qcut`. The report lists every bin's limits, count and label (e.g. `(2.5, 5]`); for k-means, the smallest and largest value, count and centre of each bin, and the share of variance between bins.

### 5.4 Recode categories

**LibreStats → Data Preparation → Coding → Recode Categories…**

Recodes or merges categories with simple rules, one per line:

```
# merge a 5-point Likert scale into three groups
1; 2 -> low
3 -> neutral
4; 5 -> high
```

```
yes; Y; y -> 1
no; N; n -> 0
* -> other
```

- `old -> new`: several old values are separated by `;`.
- `* -> value` recodes every value not listed in another rule.
- Lines starting with `#` are comments.
- Numbers match numerically (`2` matches 2 and 2.0); text matches exactly (case-sensitive, surrounding spaces ignored). New values that look like numbers are written as numbers.
- Values without a rule are **kept unchanged** (listed in the report); empty cells stay empty.
- The button **List categories in the range** fills the box with `value -> value` lines for every category found, ready to edit. It runs automatically when the dialog opens.

New column: `X recoded`.

### 5.5 Ranks and normal scores

**LibreStats → Data Preparation → Coding → Ranks and Normal Scores…**

**Ranks** (1 = smallest), new column `rank(X)`, with a choice of tie handling:

| Ties | Example for 10, 20, 20, 30 | Calc equivalent |
|---|---|---|
| average (default) | 1, 2.5, 2.5, 4 | `RANK.AVG` (ascending) |
| minimum | 1, 2, 2, 4 | `RANK.EQ` (ascending) |
| maximum | 1, 3, 3, 4 | |
| dense | 1, 2, 2, 3 | |

**Normal scores** (rank-based inverse normal transformation), new column `nscore(X)`: Φ⁻¹((r − a)/(n + 1 − 2a)), with r the average rank:

| Method | a | Formula |
|---|---|---|
| Blom (default) | 3/8 | Φ⁻¹((r − 3/8)/(n + 1/4)) |
| van der Waerden | 0 | Φ⁻¹(r/(n + 1)) |
| Tukey | 1/3 | Φ⁻¹((r − 1/3)/(n + 1/3)) |
| Rankit | 1/2 | Φ⁻¹((r − 1/2)/n) |

Normal scores give any continuous variable an approximately standard normal distribution while keeping the order of the values.

---

## 6. Data visualization

Charts are built with Calc's own chart engine, so they can be edited, copied and saved like any other Calc chart. Each tool writes the table the chart is drawn from, and the chart is linked to it: if you change a number in the table, the chart updates.

### 6.1 Histogram

**LibreStats → Data Visualization → Histogram…**

Counts how many values fall into each of a set of equal-width bins and draws them as a column chart with no gap between the bars.

All selected cells are **pooled into one list of values**, so the selection can be one column or a block of several. Text and empty cells are ignored (the number of text cells is reported). If *First row contains variable names* is ticked, the first row is skipped, and when it holds a single name, that name is used for the automatic titles.

#### Dialog

| Option | Choices | Default |
|---|---|---|
| Bins | Automatic (Sturges' rule) / Number of bins / Bin width | automatic |
| Range | Minimum and maximum; empty = the data's own minimum and maximum | empty |
| Y axis | Count / Percent | count |
| Title, X axis title, Y axis title | any text; empty = automatic ("Histogram of *name*", *name* or "Value", "Count" or "Percent") | empty |
| Bar colour | blue, light blue, green, orange, red, burgundy, grey, black | blue |
| Output | new sheet "Histogram" / to the right of the input range | new sheet |

#### Bins

With range [min, max] (the data's or your own):

| Choice | Number of bins k | Bin width |
|---|---|---|
| Automatic | Sturges (1926): k = ⌈log₂ n⌉ + 1 | (max − min) / k |
| Number of bins | k as entered | (max − min) / k |
| Bin width | k = ⌈(max − min) / width⌉ | as entered; the last bin may end beyond max |

- Bins are **half-open intervals [a, b)**: a value exactly on an edge is counted in the bin that **starts** there. The last bin also includes its upper limit, so the maximum is counted. This is the convention of `numpy.histogram`.
- **All values identical:** a single bin (with *Bin width*, from the value to value + width).
- **Your own range:** values outside [min, max] are not counted; their number is reported. Percentages are then percentages of the values counted.
- At most 1000 bins.
- Sturges' rule assumes roughly normal data and gives few bins for large samples (11 bins for n = 1000). For skewed or large data, try a fixed number of bins or a bin width.

#### Output

A two-column table (*Bin* labelled `a–b`, and *Count* or *Percent*), followed by a summary (values counted, how the bins were chosen, bin width, range, ignored text cells, values outside the range), and a column chart to the right of the table. When writing to the right of the input range, the table starts one empty column after the data; this is refused if those cells are not empty.

---

## Automation (macros)

Every analysis can be run without a dialog through the UNO service `org.librestats.Handler` (interface `com.sun.star.task.XJob`). `execute()` takes named values and returns the name of the sheet it created.

### Common arguments

| Name | Value |
|---|---|
| `Command` | see below |
| `Range` | e.g. `"$Sheet1.$A$1:$C$41"` |
| `HasLabels` | `True` (default) / `False` |
| `Document` | the spreadsheet (optional; default: the current document) |

### Commands and their arguments

| Command | Arguments (defaults) |
|---|---|
| `outliers` | `Methods` = `"tukey,mad,grubbs,dixon,mahalanobis"` (`"tukey,mad"`), `TukeyK` (1.5), `MadThreshold` (2.5), `Alpha` (0.05), `GrubbsIterate` (False), `MahalanobisAlpha` (0.001), `Highlight` (False) |
| `missing_diagnose` | `LittleTest` (True) |
| `impute` | `Method` = `mean`/`median`/`knn`/`mice` (mean), `K` (5), `Standardize` (True), `MiceMethod` = `pmm`/`norm`, `M` (5), `Iterations` (10), `Seed` (12345; empty = random), `InPlace` (False), `Highlight` (True) |
| `transform` | `method` = `log`/`sqrt`/`boxcox`/`yeojohnson`/`logit`, `base` = `e`/`10`/`2`, `offset` = `none`/`one`/`half_min`/number, `lambda` (empty = estimate), `percent`, `squeeze` (True) |
| `scale` | `method` = `centre`/`z`/`range`/`pareto` |
| `logratio` | `method` = `clr`/`ilr`, `zeros` = `error`/`replace` |
| `dummies` | `drop_first` |
| `presence` | `threshold` (0), `empty_absent` |
| `bins` | `method` = `width`/`quantile`/`kmeans`, `k` (4) |
| `recode` | `rules` (text, lines separated by newlines) |
| `ranks` | `ranks` (True), `ties` = `average`/`min`/`max`/`dense`, `scores`, `score_method` = `blom`/`vdw`/`tukey`/`rankit` |
| `histogram` | `Bins` = `auto`/`count`/`width` (auto), `K` (10), `Width`, `Min`, `Max` (empty = data), `Percent` (False), `Title`, `XTitle`, `YTitle` (empty = automatic), `Color` (RGB number, e.g. `&H004586`), `Beside` (False) |

The tools in sections 4 and 5 also accept `Beside` (True = write next to the input range).

### Example: LibreOffice Basic

```basic
Sub LogTransform
    Dim args(3) As New com.sun.star.beans.NamedValue
    args(0).Name = "Command" : args(0).Value = "transform"
    args(1).Name = "Range"   : args(1).Value = "$Sheet1.$A$1:$C$50"
    args(2).Name = "method"  : args(2).Value = "log"
    args(3).Name = "offset"  : args(3).Value = "one"
    MsgBox "Results on sheet: " & CreateUnoService("org.librestats.Handler").execute(args())
End Sub
```

### Example: Python macro

```python
from com.sun.star.beans import NamedValue

def flag_outliers(*args):
    ctx = XSCRIPTCONTEXT.getComponentContext()
    handler = ctx.ServiceManager.createInstanceWithContext("org.librestats.Handler", ctx)
    handler.execute((
        NamedValue("Command", "outliers"),
        NamedValue("Document", XSCRIPTCONTEXT.getDocument()),
        NamedValue("Range", "$Sheet1.$A$1:$C$41"),
        NamedValue("Methods", "tukey,mad,mahalanobis"),
        NamedValue("Highlight", True),
    ))
```

---

## 8. Limitations

- **Normality assumptions:** Grubbs, Dixon, Mahalanobis (χ² cutoff), Little's test and MICE's regression models assume (multivariate) normality; this is stated in the reports. Consider transforming skewed data first (section 4.1).
- **Mahalanobis distance** uses the classical mean and covariance, which outliers themselves can distort (masking). A robust version (MCD; Leys et al., 2018) is not yet available.
- **No pooling of MICE results yet:** LibreStats creates the imputed datasets but cannot yet combine analyses with Rubin's rules.
- **MAR vs MNAR** cannot be distinguished by Little's test or any other test on the observed data.
- **Performance:** pure Python is slower than compiled libraries. Data-preparation tools handle typical datasets (thousands of rows) easily; MICE with many variables, imputations and iterations takes longer.
- **Platforms:** developed and tested on Linux with LibreOffice 26.8; not yet tested on Windows or macOS.

---

## References

- Aitchison, J. (1982). The statistical analysis of compositional data. *Journal of the Royal Statistical Society B*, 44(2), 139–177.
- Blom, G. (1958). *Statistical Estimates and Transformed Beta-Variables*. Wiley.
- Box, G. E. P., & Cox, D. R. (1964). An analysis of transformations. *Journal of the Royal Statistical Society B*, 26(2), 211–252.
- Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society B*, 39(1), 1–38.
- Dixon, W. J. (1951). Ratios involving extreme values. *Annals of Mathematical Statistics*, 22(1), 68–78.
- Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal, C. (2003). Isometric logratio transformations for compositional data analysis. *Mathematical Geology*, 35(3), 279–300.
- Filzmoser, P., Hron, K., & Templ, M. (2018). *Applied Compositional Data Analysis*. Springer.
- Grubbs, F. E. (1969). Procedures for detecting outlying observations in samples. *Technometrics*, 11(1), 1–21.
- Hoaglin, D. C., & Iglewicz, B. (1987). Fine-tuning some resistant rules for outlier labeling. *Journal of the American Statistical Association*, 82(400), 1147–1149.
- Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. *Journal of Experimental Social Psychology*, 49(4), 764–766.
- Leys, C., Klein, O., Dominicy, Y., & Ley, C. (2018). Detecting multivariate outliers: Use a robust variant of the Mahalanobis distance. *Journal of Experimental Social Psychology*, 74, 150–156.
- Little, R. J. A. (1988). A test of missing completely at random for multivariate data with missing values. *Journal of the American Statistical Association*, 83(404), 1198–1202.
- Mahalanobis, P. C. (1936). On the generalised distance in statistics. *Proceedings of the National Institute of Sciences of India*, 2, 49–55.
- Martín-Fernández, J. A., Barceló-Vidal, C., & Pawlowsky-Glahn, V. (2003). Dealing with zeros and missing values in compositional data sets using nonparametric imputation. *Mathematical Geology*, 35(3), 253–278.
- Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.
- Smithson, M., & Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. *Psychological Methods*, 11(1), 54–71.
- Sturges, H. A. (1926). The choice of a class interval. *Journal of the American Statistical Association*, 21(153), 65–66.
- Tabachnick, B. G., & Fidell, L. S. (2013). *Using Multivariate Statistics* (6th ed.). Pearson.
- Troyanskaya, O., Cantor, M., Sherlock, G., Brown, P., Hastie, T., Tibshirani, R., Botstein, D., & Altman, R. B. (2001). Missing value estimation methods for DNA microarrays. *Bioinformatics*, 17(6), 520–525.
- Tukey, J. W. (1977). *Exploratory Data Analysis*. Addison-Wesley.
- van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate imputation by chained equations in R. *Journal of Statistical Software*, 45(3), 1–67.
- van den Berg, R. A., Hoefsloot, H. C. J., Westerhuis, J. A., Smilde, A. K., & van der Werf, M. J. (2006). Centering, scaling, and transformations: Improving the biological information content of metabolomics data. *BMC Genomics*, 7, 142.
- Wang, H., & Song, M. (2011). Ckmeans.1d.dp: Optimal k-means clustering in one dimension by dynamic programming. *The R Journal*, 3(2), 29–33.
- Yeo, I.-K., & Johnson, R. A. (2000). A new family of power transformations to improve normality or symmetry. *Biometrika*, 87(4), 954–959.
