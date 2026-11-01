# LibreStats Documentation

User guide for **LibreStats 0.8.0**, a statistics extension for LibreOffice Calc (26.8 or newer).

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
   - [Normal Q-Q plot](#62-normal-q-q-plot)
   - [3D surface plot](#63-3d-surface-plot)
7. [Assumption checks](#7-assumption-checks)
   - [Normality (Shapiro-Wilk)](#71-normality-shapiro-wilk)
   - [Equal variances (Levene / Brown-Forsythe)](#72-equal-variances-levene--brown-forsythe)
8. [Comparison tests](#8-comparison-tests)
   - [t-tests](#81-t-tests)
   - [ANOVA](#82-anova)
   - [Post hoc tests](#83-post-hoc-tests)
9. [Design of experiments (DOE / RSM)](#9-design-of-experiments-doe--rsm)
   - [Screening designs](#91-screening-designs)
   - [Response surface designs](#92-response-surface-designs)
   - [Analyse a screening design](#93-analyse-a-screening-design)
   - [Analyse a response surface design](#94-analyse-a-response-surface-design)
10. [Automation (macros)](#automation-macros)
11. [Limitations](#11-limitations)
12. [References](#references)

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
                                ▸ Normal Q-Q Plot…
                                ▸ 3D Surface Plot…
           ▸ Assumption Checks ▸ Normality (Shapiro-Wilk)…
                               ▸ Equal Variances (Levene / Brown-Forsythe)…
           ▸ Comparison Tests ▸ t-Tests ▸ One-Sample t-Test…
                                        ▸ Welch's t-Test (Two Samples)…
                                        ▸ Paired t-Test…
                              ▸ ANOVA ▸ One-Way ANOVA…
                                      ▸ Welch's ANOVA…
                                      ▸ Repeated-Measures ANOVA…
                              ▸ Post Hoc Tests ▸ Tukey HSD…
                                               ▸ Games-Howell…
                                               ▸ Dunnett (Each Group vs Control)…
                                               ▸ Dunn (after Kruskal-Wallis)…
           ▸ DOE / RSM ▸ New Design ▸ Screening ▸ Full Factorial (2 Levels)…
                                                ▸ Fractional Factorial 2^(k-p)…
                                                ▸ Plackett-Burman…
                                    ▸ Response Surface ▸ Central Composite (CCD)…
                                                       ▸ Box-Behnken…
                       ▸ Analyse Screening Design…
                       ▸ Analyse Response Surface Design…
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
| Normal Q-Q plot | ignored | ignored |
| Shapiro-Wilk, Levene / Brown-Forsythe | ignored | ignored, reported |
| One-sample and Welch t-test, one-way and Welch ANOVA | ignored | ignored, reported |
| Paired t-test, repeated-measures ANOVA | the whole row (subject) is left out, reported | the whole row is left out, reported |

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

The one exception is the 3D surface plot: Calc has no surface chart type, so the surface is drawn by LibreStats as a picture (SVG) and is not linked to the table.

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

### 6.2 Normal Q-Q plot

**LibreStats → Data Visualization → Normal Q-Q Plot…**

Plots the sorted values of a sample against the quantiles they would have if the data came from a normal distribution. If the data are normal, the points lie close to a straight line. A Q-Q plot shows *how* the data deviate from normality (skew, heavy tails, outliers, rounding), which a single test result such as Shapiro-Wilk (section 7.1) cannot.

#### Data layout

The same layout as the assumption checks (chapter 7):

- **One plot per column:** select one or more columns; each column is one sample. With *First row contains variable names* ticked, the names are used as titles.
- **One plot per group:** select **one column of values** and enter a **group range**: one column of group labels covering the same rows. Values with an empty group label are not used.

Text and empty cells are ignored. A sample needs at least 3 values; samples with fewer are listed under "Not plotted" and the others are still drawn.

#### Dialog

| Option | Choices | Default |
|---|---|---|
| Group labels | a range with one group label per row (optional) | empty: one plot per column |
| Plotting positions | Blom (r − 3/8) / (n + 1/4); van der Waerden r / (n + 1); Tukey (r − 1/3) / (n + 1/3); Rankit (r − 1/2) / n | Blom |
| Reference line | through the quartiles (as R's `qqline`); mean and SD; none | quartiles |
| Title, X axis title, Y axis title | any text; empty = automatic ("Normal Q-Q plot of *name*", "Theoretical quantiles (standard normal)", *name*) | empty |
| Point colour | blue, light blue, green, orange, red, burgundy, grey, black | blue |

#### Method

The i-th smallest value is plotted against the standard normal quantile Φ⁻¹((i − a) / (n + 1 − 2a)), where the plotting-position constant a is 3/8 (Blom), 0 (van der Waerden), 1/3 (Tukey) or 1/2 (Rankit). The choice hardly changes the shape of the plot; Blom is the usual default. Tied values get consecutive positions, as in R's `qqnorm`. R's `qqnorm` uses Blom's a = 3/8 for n ≤ 10 and a = 1/2 (Rankit) for larger samples; choose Rankit to match R exactly for n > 10.

| Reference line | Intercept and slope | When to use |
|---|---|---|
| Through the quartiles | the line through (z₀.₂₅, Q1) and (z₀.₇₅, Q3), with Q1 and Q3 the sample quartiles (R's default definition, type 7): slope = (Q3 − Q1) / (z₀.₇₅ − z₀.₂₅), intercept = (Q1 + Q3) / 2 | the default: it follows the middle of the data, so tails and outliers do not pull it |
| Mean and SD | intercept = mean, slope = SD (n − 1) | the normal distribution with the sample's own mean and SD; extreme values affect it |

#### Output: sheet "Q-Q plot"

For each sample: a summary (name, n, plotting positions, reference line with intercept and slope, notes), a table (*Normal quantile*, the sorted values and the *Reference line*) and an XY chart, built with Calc's chart engine and linked to the table. With several samples, the charts are stacked to the right of the tables.

#### Reading a Q-Q plot

| Pattern (sample values on the vertical axis) | Meaning |
|---|---|
| Points close to the line | consistent with a normal distribution |
| Curve bending upward at both ends (U-shape, convex) | right skew: a long tail of large values |
| Curve bending downward at both ends (∩-shape, concave) | left skew: a long tail of small values |
| Points below the line on the left and above it on the right (S-shape, steeper at the ends) | heavy tails: more extreme values than a normal distribution |
| Points above the line on the left and below it on the right | light tails (e.g. uniform-like data) |
| Horizontal steps | tied or rounded values |
| A single point far from the line at one end | a possible outlier (section 2) |

Some wiggle at the ends is normal for small samples. Judge the overall pattern, not single points.

**Checked against:** scipy's normal quantiles for all four plotting positions; R's `qnorm(ppoints(5))`; the quartile line as R's `qqline`.

### 6.3 3D surface plot

**LibreStats → Data Visualization → 3D Surface Plot…**

Draws a response z over two variables x and y as a coloured 3-D surface, e.g. extraction yield over temperature and time.

#### Data layout

Select **exactly three columns: x, y, z**, one measurement per row. With *First row contains variable names* ticked, the names are used for the axis titles.

- Rows without three numbers (empty or text cells) are skipped; their number is reported.
- Rows with the same (x, y) are **replicates**: their z values are averaged, and the table shows how many were averaged (*n*).
- x and y values that differ only by rounding noise (e.g. 0.1 × 3 and 0.3) count as the same value.

#### How the points are joined

| Points | Surface | Between the points |
|---|---|---|
| On a rectangular grid (every x with every y, as in a full factorial design; gaps allowed) | one tile per grid cell; a cell missing a corner is left open | bilinear |
| Scattered (e.g. a central composite design) | Delaunay triangulation (Bowyer, 1981; Watson, 1981) of the points, scaled to the same range | linear (flat triangles) |

The points are treated as a grid when every point is the corner of at least one complete grid cell; otherwise they are triangulated. The surface always passes **exactly through the (averaged) data**: nothing is smoothed and no model is fitted. With few scattered points (a central composite design has 9 distinct points) the surface is therefore angular. Fitting a quadratic response-surface model is a separate analysis.

#### Dialog

| Option | Choices | Default |
|---|---|---|
| Rotation | degrees, turns the surface around its vertical axis | −60 |
| Elevation | 0 (from the side) to 90 (from straight above) degrees | 30 |
| Mesh lines | draw the edges of the grid cells or triangles | on |
| Title, X, Y, Z axis title | any text; empty = automatic (column names, title "*z* by *x* and *y*") | empty |

#### Drawing

x, y and z are each scaled to [−1, 1], so a variable with large numbers does not flatten the others. The scaled box is turned by the rotation angle, tilted by the elevation and projected orthographically (no perspective). Tiles are split into small pieces (on the same surface) and drawn from the farthest to the nearest, so nearer parts hide farther ones (painter's algorithm; Newell et al., 1972). Each piece is coloured by its height on the **viridis** scale (van der Walt & Smith, 2015), which is perceptually uniform and readable with colour-vision deficiencies and in greyscale. Axis ticks and the colour legend are in the original units.

#### Output

A new sheet "3D Surface" with a summary (points plotted, grid or triangulation, view angles, skipped rows, notes), the table of plotted points (x, y, mean z, n), and the picture to the right.

- The picture is **not linked** to the data: change the data or the angles and run the tool again. To see another side, change the rotation (e.g. −60, 30, 120, 210).
- It is an ordinary Calc image: move, resize, copy it, or right-click → *Save…* to export it as SVG or PNG.
- Limits: a triangulation of at most 2000 scattered points, at most 20 000 tiles.

---

## 7. Assumption checks

Tests of the two assumptions behind the classical t-tests and ANOVA (chapter 8): that the data in each group are **normally distributed**, and that all groups have the **same variance**. Each report gives the test per variable or group, a conclusion at the chosen significance level (0.05, 0.01 or 0.10), notes and references.

**Data layout** (the same as for the independent-samples tests in chapter 8):

- **Samples as columns:** each **column** of the input range is one variable (normality) or one group (equal variances). With *First row contains variable names* ticked, the first row gives the names; otherwise columns are called "Column B", etc.
- **Values and group labels:** the input range is **one column of values**, and the **group range** is one column of group labels covering **the same rows**. Labels can be text or numbers; groups appear in the order of their first appearance. Values with an empty group label are not used (their number is reported).

Text cells in the values are ignored and counted in the report; empty cells are ignored.

### 7.1 Normality (Shapiro-Wilk)

**LibreStats → Assumption Checks → Normality (Shapiro-Wilk)…**

Tests the null hypothesis that the data come from a normal distribution, for each column or group separately.

#### Dialog

| Option | Meaning | Default |
|---|---|---|
| Input range | the data (columns, or one column of values) | current selection |
| Group labels | optional; one group label per row of the input range | empty: test each column |
| Significance level | 0.05, 0.01 or 0.10 | 0.05 |

#### Method

The Shapiro-Wilk statistic (Shapiro & Wilk, 1965) is

W = (Σ aᵢ x₍ᵢ₎)² / Σ (xᵢ − x̄)²

with x₍ᵢ₎ the sorted values and aᵢ coefficients derived from the expected normal order statistics. W is at most 1; values well below 1 indicate non-normality. The coefficients and the p-value use Royston's approximation, algorithm AS R94 (Royston, 1992, 1995). This is the same method as R's `shapiro.test` and scipy's `shapiro`. For n = 3 the p-value is exact.

- **3 to 5000 values** per variable or group. Larger samples are still tested, with a note that the p-value is approximate beyond n = 5000.
- **Not computed** (reported as such in the table, the other variables are still tested): fewer than 3 values, or all values equal.
- **Tied values** are noted: heavily rounded data (few distinct values) can make the test reject normality although the underlying quantity is normal.

#### Output: sheet "Normality"

The input range, the null hypothesis and the significance level, then a table with one row per variable or group: **n, mean, SD, W, p-value** and the **conclusion**, either "normality rejected (p < α)" or "no evidence against normality". Then notes and references.

#### Interpretation

- A non-significant result does **not prove** normality; it means the data do not contradict it. With small samples the test has little power and misses real deviations.
- With large samples it detects deviations far too small to matter for a t-test or ANOVA.
- Always look at the data too: a **Normal Q-Q plot** (section 6.2) or a histogram (section 6.1) shows what kind of deviation there is. Strong skew and outliers matter more for tests of means than a small but significant departure.
- If normality is clearly violated, a transformation (section 4.1, e.g. log or Box-Cox) often helps. For comparing groups, Dunn's test after Kruskal-Wallis (section 8.3) does not assume normality.

**Checked against:** scipy's `shapiro` for n = 3 to 5000 on normal, skewed, uniform and heavy-tailed data (W and p within 10⁻¹⁰), and the exact result for n = 3.

### 7.2 Equal variances (Levene / Brown-Forsythe)

**LibreStats → Assumption Checks → Equal Variances (Levene / Brown-Forsythe)…**

Tests the null hypothesis that all groups have the same variance. The groups are the columns of the input range, or the groups given by the group range.

#### Dialog

| Option | Meaning | Default |
|---|---|---|
| Input range | the data (one column per group, or one column of values) | current selection |
| Group labels | optional; one group label per row of the input range | empty: each column is a group |
| Significance level | 0.05, 0.01 or 0.10 | 0.05 |
| Deviations from the group | **median**: the Brown-Forsythe version (recommended); **mean**: the original Levene test | median |

#### Method

Each value is replaced by its absolute deviation from its group's centre, zᵢⱼ = |xᵢⱼ − centreᵢ|, and the zᵢⱼ are compared by a one-way ANOVA:

F = [(N − k) / (k − 1)] · Σ nᵢ (z̄ᵢ − z̄)² / Σᵢ Σⱼ (zᵢⱼ − z̄ᵢ)²,  with k − 1 and N − k degrees of freedom

(k groups, N values in total). If the groups have different spreads, the mean absolute deviations differ and F is large.

- **Centre = median** (Brown & Forsythe, 1974): keeps its error rate when the data are skewed or heavy-tailed. This is the default here, as in R's `car::leveneTest`.
- **Centre = mean** (Levene, 1960): the original test; accurate only for normal data.

Requirements: at least **2 groups**, each with at least **2 values**. If the absolute deviations do not vary within any group (e.g. every group has exactly two values), the test cannot be computed and a message says so.

#### Output: sheet "Equal variances"

The input range, the null hypothesis, the centre and the significance level, then:

- a table per group: **n, mean, median, SD, variance**;
- the test: **F, df1, df2, p-value** and the **conclusion**, either "variances differ (p < α)" or "no evidence of unequal variances";
- the **ratio of the largest to the smallest variance**, a direct measure of how different the spreads are;
- notes and references.

#### Interpretation

- As with any test, a non-significant result does not prove equal variances, and with large groups small differences become significant. The variance ratio shows whether a difference is large.
- Unequal variances matter most when group sizes differ. The tests in chapter 8 that do **not** assume equal variances are Welch's t-test, Welch's ANOVA and the Games-Howell post hoc test. They are good default choices whatever this test says. The classical one-way ANOVA and Tukey's HSD do assume equal variances.

**Checked against:** scipy's `levene` with `center="median"` and `center="mean"` (F and p within 10⁻¹⁰), and a hand-computed example.

---

## 8. Comparison tests

Tests of whether means differ. Every report gives the descriptive statistics (n, mean, SD, SE), the test, an effect size, a conclusion at the chosen significance level (0.05, 0.01 or 0.10), notes on the assumptions, and references.

There are two data layouts:

- **Independent samples** (one-sample and Welch t-test, one-way and Welch ANOVA): each **column** of the input range is one sample, or, if a **group range** is given, the input range is one column of values and the group range (same rows) holds the group labels. This is the same layout as the assumption checks (chapter 7).
- **Related samples** (paired t-test, repeated-measures ANOVA): **one column per condition and one row per subject**. Only complete rows are used; rows with an empty or text cell in any selected column are left out, and their number is reported.

Normality can be checked with **Assumption Checks → Normality**, and equal variances with **Assumption Checks → Equal Variances** (chapter 7).

### 8.1 t-tests

**LibreStats → Comparison Tests → t-Tests**

| Test | Data | Statistic | df | Effect size |
|---|---|---|---|---|
| One-sample (Student, 1908) | one or more samples, each tested against μ0 | t = (x̄ − μ0) / (s / √n) | n − 1 | d = (x̄ − μ0) / s |
| Welch (1947) | exactly 2 samples (2 columns, or 2 group labels) | t = (x̄₁ − x̄₂) / √(s₁²/n₁ + s₂²/n₂) | Welch-Satterthwaite | d = (x̄₁ − x̄₂) / √((s₁² + s₂²) / 2) (Cohen, 1988) |
| Paired | exactly 2 columns, one row per pair | t = d̄ / (s_d / √n) on the differences d = first − second | n − 1 | d_z = d̄ / s_d (Lakens, 2013) |

- **Alternative hypothesis:** two-sided (≠), or one-sided greater (>) or less (<). Differences are always *first − second* (first column, or the first group label in the data).
- **Confidence interval** at 1 − α: for the mean (one-sample), the difference of means (Welch) or the mean difference (paired). For a one-sided test it is one-sided, with an infinite bound (shown as ∞ or −∞), as in R.
- The paired report also shows the correlation of the pairs.
- **Why Welch only:** Welch's test keeps its error rate when the variances differ and loses almost nothing when they are equal, so it is recommended as the default two-sample test (Delacre et al., 2017). Student's pooled-variance test is therefore not offered.

### 8.2 ANOVA

**LibreStats → Comparison Tests → ANOVA**

| Test | Data | Output |
|---|---|---|
| One-way ANOVA (Fisher, 1925) | 2 or more independent groups | ANOVA table (SS, df, MS, F, p) for between groups, within groups and total; η² and ω² (Olejnik & Algina, 2003) |
| Welch's ANOVA (Welch, 1951) | 2 or more independent groups, each with ≥ 2 values that vary | F, df1 = k − 1, df2 (not an integer), p |
| Repeated-measures ANOVA | 2 or more conditions (columns), one row per subject | ANOVA table for conditions, subjects, error and total; Mauchly's test; corrected p-values; partial η² and generalized η² (Bakeman, 2005) |

- The classical one-way ANOVA assumes **equal variances**; when they differ, especially with unequal group sizes, use Welch's ANOVA (Delacre et al., 2019). With 2 groups Welch's ANOVA equals Welch's t-test (F = t²).
- **Repeated measures: sphericity.** For 3 or more conditions the report gives **Mauchly's (1940) test** on orthonormal contrasts (it needs more subjects than conditions − 1) and the F-test with four sets of degrees of freedom: sphericity assumed, **Greenhouse-Geisser** ε (1959), **Huynh-Feldt** ε (1976; with one group of subjects the Lecoutre (1991) correction is identical; capped at 1) and the lower bound 1/(k − 1). Each corrected test uses ε × df1 and ε × df2. Mauchly's test has little power in small samples, so many authors report a corrected p-value regardless; a common rule is Huynh-Feldt when the Greenhouse-Geisser ε > 0.75, otherwise Greenhouse-Geisser (Girden, 1992).
- With 2 conditions sphericity always holds and the repeated-measures F equals the squared paired t.

### 8.3 Post hoc tests

**LibreStats → Comparison Tests → Post Hoc Tests**

Multiple comparisons between groups, with the same data layout as the one-way ANOVA (columns, or values plus a group range). p-values and confidence intervals are **adjusted** for the whole set of comparisons, so the chance of at least one false positive stays at α (Dunn's test with Benjamini-Hochberg: the false discovery rate).

| Test | Compares | Statistic | Distribution | Assumes |
|---|---|---|---|---|
| Tukey HSD (Tukey, 1953; Tukey-Kramer for unequal n, Kramer, 1956) | all pairs of means | q = \|x̄ᵢ − x̄ⱼ\| / √(MSE/2 · (1/nᵢ + 1/nⱼ)), MSE pooled over all groups, df = N − k | studentized range (k, df) | normal data, equal variances |
| Games-Howell (Games & Howell, 1976) | all pairs of means | t = (x̄ᵢ − x̄ⱼ) / √(sᵢ²/nᵢ + sⱼ²/nⱼ), Welch df per pair | studentized range (k, Welch df), q = √2·\|t\| | normal data; variances may differ |
| Dunnett (1955) | each group with one control | t = (x̄ᵢ − x̄₀) / √(MSE (1/nᵢ + 1/n₀)), df = N − k | joint multivariate t of the k − 1 comparisons (correlations λᵢλⱼ, λᵢ = √(nᵢ/(nᵢ + n₀))) | normal data, equal variances |
| Dunn (1964) | all pairs of mean ranks | z = (R̄ᵢ − R̄ⱼ) / √(σ² (1/nᵢ + 1/nⱼ)), σ² = N(N+1)/12 − Σ(t³ − t)/(12(N − 1)) | standard normal, then adjusted | none on the distribution shape |

- **Tukey HSD** and **Games-Howell** report every pair with the difference (*Group 1 − Group 2*), its standard error, the statistic, the adjusted p-value and the simultaneous confidence interval. Tukey's HSD does not need a significant ANOVA first. Use Games-Howell when the variances differ (it can be slightly liberal with fewer than about 6 values per group).
- **Dunnett:** enter the **control group** as a column name or group label (empty = the first group). The test can be two-sided or one-sided (each group greater, or less, than the control); one-sided confidence intervals have an infinite bound. It is more powerful than Tukey when only the comparisons with the control matter.
- **Dunn:** the report first gives the **Kruskal-Wallis** test (tie-corrected H, df = k − 1, χ² p-value) and the median and mean rank of each group, then all pairs with unadjusted and adjusted p-values. Adjustments: **Holm** (default; step-down, controls the family-wise error rate, uniformly more powerful than Bonferroni), **Bonferroni** (Dunn's original), **Šidák**, or **Benjamini-Hochberg** (false discovery rate).
- **Computation:** the studentized range and Dunnett distributions are computed by numerical integration (Gauss-Legendre quadrature over the scale of the variance estimate and over the normal), giving the upper tail directly so that small p-values stay accurate; critical values are found by root finding. This takes a fraction of a second per comparison, so Games-Howell with many groups (e.g. 3 s for 8 groups) is the slowest tool.

---

## 9. Design of experiments (DOE / RSM)

**LibreStats → DOE / RSM → New Design → Screening → …** and **→ Response Surface → …**

Creates a **design sheet**: the list of runs to carry out, in random order, with the factor settings in coded and actual units and empty columns for the responses you measure.

### 9.1 Screening designs

Screening designs test many factors at two levels each (low and high) to find the few that matter.

| Menu item | Design | Runs |
|---|---|---|
| Full Factorial (2 Levels)… | every combination of low and high levels (2 to 8 factors) | 2^k |
| Fractional Factorial 2^(k-p)… | a regular fraction of the full factorial (3 to 11 factors) | 4 to 128 |
| Plackett-Burman… | Plackett & Burman (1946) designs (2 to 23 factors) | 8, 12, 16, 20 or 24 |

The three menu items open the same dialog with the design type preselected.

#### Factor table

Before opening the dialog, type the factors on a sheet, one per row, and select the table:

| Name | Unit | Low | High |
|---|---|---|---|
| Temperature | °C | 40 | 80 |
| Time | min | 10 | 30 |
| Catalyst | | Pd | Pt |

- *Unit* may be empty, or the column left out (Name, Low, High).
- Levels can be **text** (a categorical factor such as a catalyst type). Centre points then are not possible.
- Factors get the letters A, B, C, … in table order (I is skipped, as usual, because it stands for the identity in alias lists).
- Low is coded −1 and High +1; centre points (coded 0) are run at the midpoint (Low + High) / 2. Section 8.2 explains the conversion in detail.

#### Dialog

| Option | Meaning | Default |
|---|---|---|
| Design | full factorial / fractional factorial / Plackett-Burman | from the menu item |
| Runs | fractional: the size of the fraction; Plackett-Burman: 8, 12, 16, 20 or 24. Empty: automatic (see below) | automatic |
| Centre points per block | runs with every numeric factor at its midpoint (coded 0); they check for curvature and estimate pure error | 0 |
| Replicates | the whole design is run this many times | 1 |
| Blocks | total number of blocks: replicates × 1, 2, 4 or 8 (see below) | 1 |
| Randomise the run order | shuffles the runs within each block | on |
| Seed | whole number for a reproducible order; empty = a new random seed, which is reported | empty |
| Response columns | names of the empty response columns, separated by commas | Response |

#### Fractional factorials

LibreStats uses the **minimum-aberration** fraction (Fries & Hunter, 1980) for each number of factors and runs: of all fractions of that size, the one with the fewest short words in its defining relation, so the fewest low-order interactions aliased with each other. The generators were found by an exhaustive search (`tools/make_ff_table.py`). Their resolutions agree with the NIST/SEMATECH e-Handbook table, and their word-length patterns with Montgomery's Table 8.14.

| Resolution (Box & Hunter, 1961) | Meaning |
|---|---|
| III | main effects are aliased with two-factor interactions |
| IV | main effects are clear of two-factor interactions; two-factor interactions are aliased with each other |
| V | main effects and two-factor interactions are clear of each other |

**Automatic size:** the smallest fraction of resolution IV or higher; if there is none (3 factors), the 4-run resolution III fraction. If you enter a size that does not exist, the message lists the sizes available with their resolutions.

The summary lists the **design generators** (e.g. E = ABC), the **defining relation** and the **alias chains** of all main effects and two-factor interactions, showing interactions up to the third order (e.g. `A + BCE + BDF`).

#### Plackett-Burman designs

Built from the cyclic first rows of Plackett & Burman (1946) plus a row of all low levels; with k factors the first k columns are used, and the unused columns can serve to estimate error. Automatic size: the smallest N ≥ k + 1.

- With 12, 20 or 24 runs the aliasing is **partial**: each two-factor interaction is partly aliased with every main effect it does not involve. Use these designs only when interactions are expected to be small.
- With 8 or 16 runs a Plackett-Burman design is a regular resolution III fraction.

#### Blocks

Blocks group runs done under similar conditions (same day, same batch of raw material). *Blocks* is the total number:

- **Blocks = replicates:** each replicate is one block.
- **Blocks = replicates × 2, 4 or 8** (factorial designs only): each replicate is split by **confounding** interactions with blocks. LibreStats chooses the block generators so that no main effect is confounded with blocks and as few low-order interactions as possible are. For full factorials this reproduces the quality of Montgomery's Table 7.8. The confounded effects (with their aliases) are listed.
- Plackett-Burman designs can only be blocked by replicate.

Centre points are added to each block. Runs are randomised within blocks; the blocks stay in order.

#### Output: sheet "Design"

| Column | Content |
|---|---|
| Std order | the run's place in the unrandomised design (Yates order, A changing fastest; replicate by replicate, block by block, centre points last in each block) |
| Run order | the order in which to carry out the runs; the sheet is sorted by it |
| Point type | 1 = factorial point, 0 = centre point (response surface designs: −1 = axial point, 2 = Box-Behnken edge midpoint) |
| Block | 1, 2, … |
| *factor* (coded) | −1, 0, +1 |
| *factor* (*unit*) | the actual setting |
| response columns | empty and shaded, to be filled in |

To the right: a summary (design, number of runs, replicates, blocks, centre points, randomisation seed), the factor table with letters, and for fractional factorials the generators, defining relation and alias chains, the effects confounded with blocks, notes and references. The header row is frozen.

### 9.2 Response surface designs

Response surface designs have three or more levels per factor, so a **full quadratic model** (linear terms, squares and two-factor interactions) can be fitted to find an optimum. All factors must be numeric.

| Menu item | Design | Factors | Runs (default centre points) |
|---|---|---|---|
| Central Composite (CCD)… | factorial (cube) points + 2k axial points at ±α on each axis + centre points (Box & Wilson, 1951) | 2 to 7 | 13, 20, 31, 32, 53, 92 for k = 2 … 7 |
| Box-Behnken… | each pair (k = 3–5) or triple (k = 6, 7) of factors at ±1 with the others at 0, + centre points (Box & Behnken, 1960); 3 levels, no corner points | 3 to 7 | 15, 27, 46, 54, 62 |

#### Factor table

Before opening the dialog, type the factors on a sheet, one per row, and select the table (or enter its range in the dialog's *Factor table* field):

| Name | Unit | Low | High |
|---|---|---|---|
| Temperature | °C | 40 | 80 |
| Time | min | 10 | 30 |
| pH | | 5 | 7 |

| Column | Required? | Content |
|---|---|---|
| Name | yes | the factor name; names must differ. Used in the column headers of the Design sheet |
| Unit | no (may be empty, or the column left out: Name, Low, High) | shown in the headers, e.g. `Temperature (°C)` |
| Low | yes | a **number**: the actual setting coded −1 |
| High | yes | a **number**: the actual setting coded +1; must differ from Low |

- Tick *First row contains column headings* if the table has a header row (as above).
- **Only numeric levels.** A factor with text levels (e.g. a catalyst type Pd/Pt) cannot be used in a response surface design; the dialog stops with a message naming it. Put such factors in a screening design (section 9.1), or make a separate response surface design for each category.
- 2 to 7 factors for a CCD, 3 to 7 for Box-Behnken. Factors get the letters A, B, C, … in table order (I is skipped).

**From actual to coded values.** Each factor is coded around the midpoint of Low and High:

coded = (actual − centre) / half-range, with centre = (Low + High) / 2 and half-range = (High − Low) / 2

so Low → −1, the midpoint → 0 and High → +1. LibreStats builds the design in coded units, then converts each run back to actual settings: actual = centre + coded × half-range. The Design sheet shows both. What Low and High mean depends on the design:

| Design | Low and High are | Example: Temperature 40–80 °C, Time 10–30 min, CCD with 2 factors (α = 1.4142) |
|---|---|---|
| CCD, rotatable / orthogonal / custom α (circumscribed, the default) | the **cube** (factorial) levels ±1. The axial points at ±α lie **outside** Low–High | Temperature 31.72, 40, 60, 80, 88.28 °C; Time 5.86, 10, 20, 30, 34.14 min (at −α, −1, 0, +1, +α) |
| CCD, *Inscribed* ticked | the **outer limits**: ±α is put at Low and High, the cube shrinks to ±1/α | Temperature 40, 45.86, 60, 74.14, 80 °C; Time 10, 12.93, 20, 27.07, 30 min |
| CCD, face-centred (α = 1) | the cube levels, which are also the extremes: three levels per factor | Temperature 40, 60, 80 °C; Time 10, 20, 30 min |
| Box-Behnken | the extremes: three levels −1, 0, +1 | Temperature 40, 60, 80 °C; Time 10, 20, 30 min; pH 5, 6, 7 |

So for a standard CCD, **choose Low and High as the cube levels, not as the limits of what is possible.** If the axial settings beyond them are impossible (e.g. a negative time, or a temperature the equipment cannot reach), tick *Inscribed* or choose the face-centred α. LibreStats warns when α > 1 puts the axial points outside Low–High, and names any factor whose setting at −α would become negative.

The actual settings on the Design sheet are rounded to about 4 significant digits of the half-range (e.g. 88.28 °C, 948.8 rpm); the coded values keep full precision. The summary to the right of the design lists each factor's settings at −α, −1, 0, +1 and +α (Box-Behnken: −1, 0, +1), so you can check them before running the experiment.

#### Dialog

The dialog has the shared options of section 9.1 (replicates, blocks, randomisation, seed, response columns) plus:

| Option | Meaning | Default |
|---|---|---|
| Axial distance α (CCD) | **rotatable**: α = n_F^(1/4) (n_F = number of cube points), so the prediction variance depends only on the distance from the centre; **orthogonal blocking**: the α that makes the blocks orthogonal (needs blocks); **face-centred**: α = 1, three levels, axial points on the faces of the cube; **custom**: any positive value | rotatable |
| Resolution V fraction as factorial part (CCD, 5 or more factors) | the cube is the smallest fraction of resolution V or higher (2^(5−1), 2^(6−1), 2^(7−1)) instead of the full 2^k; all quadratic-model terms stay estimable | on |
| Inscribed (CCD) | the axial points are put at Low and High and the cube is shrunk to ±1/α in actual units. Use it when Low and High are limits that cannot be exceeded | off |
| Centre points | empty = automatic, see below. Unblocked: per replicate; blocked: per block (CCD: per cube block) | automatic |
| Centre points in the axial block (CCD in blocks) | empty = automatic | automatic |

#### Centre points (automatic)

- **CCD, unblocked:** the uniform-precision number of Box & Hunter (1957), at which the prediction variance at the centre is about the same as at distance 1: N = λ₄ (n_F + 2α²)² / n_F with their tabulated λ₄. This gives 5, 6, 7, 6 (half fraction) / 10 (full), 9 / 15, 14 / 21 centre points for k = 2 … 7, as in the published table.
- **Box-Behnken, unblocked:** Box & Behnken's (1960) numbers: 3, 3, 6, 6, 6 for k = 3 … 7.
- **In blocks:** 2 per block (CCD: 2 per cube block and 2 in the axial block), so each block has its own repeated points.

#### Blocks

| Design | Blocks per replicate |
|---|---|
| CCD | 1; 2 = cube block + axial block; 3 = cube split in 2 + axial block; 5 = cube split in 4 + axial block. The cube is split by confounding high-order interactions with blocks; splits that would confound a two-factor interaction are refused |
| Box-Behnken | 1; for 4 factors also 3 orthogonal blocks (pairs AB + CD, AD + BC, AC + BD) |

As for screening designs, *Blocks* is the total: replicates × the blocks per replicate.

For a CCD in blocks, **orthogonal blocking** (the block effects do not bias the model coefficients) holds when n_F / (n_F + n_cf) = 2α² / (2k + n_ca), with n_cf centre points in the cube blocks and n_ca in the axial block (Box & Hunter, 1957). The *orthogonal blocking* α solves this. Examples: 2 factors in 2 blocks give α = √2, both rotatable and orthogonal; 3 factors in 3 blocks with 2 centre points each give α = 1.633; 4 factors in 3 blocks give α = 2, again both.

#### Design properties in the summary

The summary checks the design that was built rather than relying on rules of thumb:

| Property | How it is determined |
|---|---|
| Rotatable | moment conditions up to order 4: odd moments vanish, all [ii], [iiii] and [iijj] are equal, and [iiii] = 3 [iijj] (Box & Hunter, 1957). Box-Behnken designs are rotatable for 4 and 7 factors only |
| Blocks orthogonal to the quadratic model | every model column has the same mean in every block |
| Terms in the full quadratic model | 1 + 2k + k(k − 1)/2 |
| Residual degrees of freedom | runs − terms − (blocks − 1) |
| Pure-error degrees of freedom | from repeated points (centre points, replicates) in the same block; with 0, lack of fit cannot be tested |

The factor table shows the actual settings at −α, −1, 0, +1 and +α (Box-Behnken: −1, 0, +1). Actual settings are rounded to about 4 significant digits of the half-range (e.g. 948.8 rpm, 114.88 °C); coded values are kept at full precision and shown with up to 4 decimals. If α > 1, a note warns that the axial points lie outside Low–High, and names any factor whose setting at −α would be negative.

### 9.3 Analyse a screening design

**LibreStats → DOE / RSM → Analyse Screening Design…**

After the runs are done and the responses are typed into the shaded columns of a Design sheet (section 9.1), this tool estimates which factors matter.

**Only for screening designs** (full or fractional factorial, Plackett-Burman). A **response surface design** (central composite or Box-Behnken, section 9.2) has axial or edge points and needs a quadratic model: use *Analyse Response Surface Design* (section 9.4). Opening this tool on such a Design sheet shows a message saying so.

#### Dialog

| Option | Meaning | Default |
|---|---|---|
| Design table | the table on a Design sheet, header row included, laid out as described under *Required table layout* below; filled in automatically when a Design sheet is active | active Design sheet |
| Response | one of the response columns | the first |
| Model | main effects and two-factor interactions / main effects only | with interactions |
| Significance level | α for the tests, the critical t and Lenth's margins | 0.05 |
| Pareto chart and half-normal plot | draw the charts | on |

#### Required table layout

A Design sheet made with *New Design* (section 9.1) already has this layout. Fill in the shaded response columns and you are done. To analyse a design made elsewhere, build a table with the same layout, for example a 2² design with one centre point:

| Std order | Run order | Point type | Block | Temp (coded) | Time (coded) | Temp (°C) | Time (min) | Yield | Purity |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 3 | 1 | 1 | -1 | -1 | 40 | 10 | 61.2 | 97.1 |
| 2 | 1 | 1 | 1 | 1 | -1 | 80 | 10 | 70.5 | 95.4 |
| 3 | 5 | 1 | 1 | -1 | 1 | 40 | 30 | 64.8 | |
| 4 | 2 | 1 | 1 | 1 | 1 | 80 | 30 | 78.0 | 94.2 |
| 5 | 4 | 0 | 1 | 0 | 0 | 60 | 20 | 69.9 | 96.0 |

| Column | Required? | Header | Content |
|---|---|---|---|
| Std order, Run order | no (not used by the analysis) | any | anything. When cell A1 of the active sheet is `Std order`, the dialog finds the table by itself; otherwise enter its range |
| Point type | **yes** | exactly `Point type` | `1` = factorial run, `0` = centre point. Other values (`-1` axial point, `2` edge midpoint) mean a response surface design and are refused |
| Block | **yes** | exactly `Block` | a number per run; all `1` when the design is not blocked |
| coded factor columns | **yes, at least 2, side by side** | the factor name followed by ` (coded)`, e.g. `Temp (coded)` | `-1` or `1` in factorial runs, `0` in centre points; nothing else is accepted. The columns become factors A, B, C, … in this order (I is skipped) |
| actual-level columns | **yes, exactly one per factor**, directly after the last coded column | any | ignored by the analysis. They must be there because the responses are looked for after them |
| response columns | **yes, at least one**, after the actual-level columns | the response name (chosen in the dialog) | numbers. An **empty** cell means that run has no result and is left out; text is an error |

- The **first row of the range is the header row**. The table ends at the first column whose header is empty, so keep an empty column between the table and anything to its right (the design summary on a Design sheet is placed that way).
- *Point type* and *Block* may stand anywhere in the table; LibreStats puts them before the factors.
- The **rows can be in any order** (standard or run order). Completely empty rows are skipped.
- There must be **at least as many runs with a response as model terms**: the constant, one per extra block, one for centre points (if any) and one per effect. Example: the table above has 5 runs and the default model needs 5 terms (constant, centre points, A, B, AB). If one response is missing, choose *Main effects only*. With exactly as many runs as terms the model is saturated, and Lenth's method supplies the error (see below).

#### Model and effects

Effects are estimated by **least squares** on the coded −1/+1 columns. The model holds:
- a constant;
- the block terms;
- a centre-point term, when there are centre points, so that curvature does not bias the effects or the error;
- the main effects, and by default the two-factor interactions.

An **effect** is twice the coefficient: the change in the mean response when the factor goes from low to high. The table gives effect, coefficient, standard error of the effect, t, p and "significant" at α.

**Aliasing is read from the design itself:** two effect columns that are equal (or opposite) over the factorial runs are aliased, and only the first (lowest-order) effect of each alias chain is estimated. The report lists every chain, e.g. `AB + CE + DF (+ 5 higher-order)`. An effect therefore estimates the sum of its chain. Effects confounded with blocks are reported and left out. In **non-regular designs** (Plackett-Burman with 12, 20 or 24 runs), two-factor interactions are only partially aliased with main effects; they are left out of the model and a note says so.

If a response is missing, least squares still gives the effects, but they are then slightly correlated (no longer independent).

#### Error estimate

| Situation | t and p use |
|---|---|
| Residual degrees of freedom (replicates, centre points, or a model smaller than the design) | the residual mean square |
| Saturated model (no residual df), every factorial run has a response | Lenth's pseudo standard error: t = effect / PSE with m/3 df |
| Saturated and responses missing | no t or p values (a note says so) |

When there are repeated runs, the residual is split into **pure error** and **lack of fit**, with an F test.

#### Lenth's method

Lenth (1989) estimates the noise from the effects themselves, assuming most of them are inactive. It uses all m estimable contrasts of the factorial runs, e.g. all 15 of an unreplicated 2^4:

- s0 = 1.5 × median |c|, and PSE = 1.5 × median of the |c| that are smaller than 2.5 s0
- margin of error ME = t(1 − α/2; m/3) × PSE
- simultaneous margin of error SME = t(γ; m/3) × PSE, with γ = (1 + (1 − α)^(1/m)) / 2

Effects beyond ME are likely active; effects beyond SME remain so even after allowing for testing all m effects at once. Lenth's method needs a response for every factorial run. For Montgomery's Example 6.2 (filtration rate) it gives PSE = 2.625, ME = 6.75 and SME = 13.70, as in the book.

#### Curvature check

With centre points, the centre-point term estimates mean(centre) − mean(factorial), adjusted for blocks. Its sum of squares (without blocks, n_F n_C (ȳ_F − ȳ_C)² / (n_F + n_C)) is tested with F against **pure error**, which needs at least 2 centre points or replicated runs (Montgomery, 2017, section 6.8). Significant curvature means that a first-order model with interactions does not describe the region: augment the design to a response surface design (section 9.2).

#### Charts (native Calc charts)

- **Pareto chart:** horizontal bars of |t| (or Lenth's pseudo-t), largest on top. It has a red reference line at the critical t; with Lenth's method, red ME and orange SME lines.
- **Half-normal plot** (Daniel, 1959): |effect| of all contrasts against the half-normal quantiles Φ⁻¹(0.5 + 0.5 (i − 0.5)/m), with the line |effect| = PSE × quantile. Inactive effects lie near the line; effects beyond ME are labelled.

The tables the charts are drawn from are on the report sheet, to the left of the charts.

### 9.4 Analyse a response surface design

**LibreStats → DOE / RSM → Analyse Response Surface Design…**

Fits a quadratic model to the results of a response surface design (central composite or Box-Behnken, section 9.2), and reports where the response is highest or lowest. It also works on two-level factorials with centre points, but only with the model without square terms (see below).

#### Required table layout

The same layout as for the screening analysis (section 9.3, *Required table layout*): a Design sheet made with *New Design*, with the responses typed into the shaded columns. Two differences:

- The **coded values may be any numbers**, e.g. axial points at ±1.414 or ±1.682. *Point type* may be `1` (cube), `-1` (axial), `0` (centre) or `2` (Box-Behnken edge midpoint); the analysis does not use it.
- The **actual-level columns are used:** they must be **numbers** (no text levels), and must agree with the coded columns (actual = centre + half-range × coded). LibreStats reads each factor's centre and half-range from them, to report the results in actual units. The settings on a Design sheet are rounded (e.g. 31.72 °C for −α), so LibreStats recognises the rows that hold the exact Low and High you typed: the cube points in a standard CCD, or the axial points in an inscribed one. A table made by hand with unrounded values also works.

As before, the rows can be in any order, an empty response cell leaves that run out, and the design needs more runs with a response than model terms (6 for 2 factors, 10 for 3, 15 for 4).

#### Dialog

| Option | Meaning | Default |
|---|---|---|
| Design table | the table on a Design sheet (header row included); filled in automatically when a Design sheet is active | active Design sheet |
| Response | one of the response columns | the first |
| Model | **full quadratic** (linear, squares, two-factor interactions) / linear + two-factor interactions / linear | full quadratic |
| Significance level | α for the tests | 0.05 |
| Response map and 3D surface | draw them, for the chosen *horizontal* and *vertical* factor; all other factors are held at their centre | on; first two factors |

#### Model

In coded units: y = b0 + Σ bᵢxᵢ + Σ bᵢᵢxᵢ² + Σ bᵢⱼxᵢxⱼ, fitted by least squares. Blocks are included with sum-to-zero coding, so the constant and every prediction refer to an **average block**.

A two-level design with centre points has only three levels where all factors are at −1, 0 or +1 together. It cannot separate the square terms, and the full quadratic model is refused with a message. Use *linear + two-factor interactions* there, or add axial points.

#### Report

| Part | Content |
|---|---|
| Fit | S (residual standard deviation), R², adjusted R², **predicted R²** (from PRESS: how well the model predicts runs left out one at a time) |
| Coefficients (coded units) | coefficient, SE, t, p and "significant" for every term. Coded coefficients can be compared with each other: a larger \|b\| means a larger effect over the design range |
| Analysis of variance | Model, Linear, Square, 2-way interaction and Blocks with adjusted sums of squares (the increase in residual SS when the group is removed), F and p; residual split into **lack of fit** and **pure error** (repeated points, mostly centre points) with an F test. Significant lack of fit means the quadratic model does not describe the response well enough |
| Equation in actual units | the same model with the factors in their own units (for an average block) |
| Factor levels | each factor's centre and half-range, and the coded range of the design |
| Stationary point and canonical analysis | the point where the slope is zero, x_s = −B⁻¹b / 2, in coded and actual units, with its predicted response; the **eigenvalues** of B: all negative = **maximum**, all positive = **minimum**, mixed signs = **saddle point**. An eigenvalue near 0 means a **ridge**: settings along its direction give nearly the same response. The report says whether the point lies inside the experimental region |
| Maximum / minimum within the experimental region | the exact best and worst predicted response within the box of the coded settings that were run, with the settings. For a maximum inside the region this is the stationary point; for a saddle point or a stationary point outside the region it is on the boundary |

The **experimental region** is the box from the lowest to the highest coded setting of each factor (for a CCD, ±α). For a rotatable CCD the corners of this box lie beyond the axial points, so predictions there are extrapolations. Prefer optima that lie within the region actually explored, and confirm any optimum with new runs.

#### Response map and 3D surface

- **Response map:** a 21 × 21 grid of predicted values over the two chosen factors, in actual units: the horizontal factor across, the vertical factor down, high values at the top. Each cell is coloured by its value (viridis: dark blue = lowest, yellow = highest on the map). It is ordinary cells, so you can read and copy the numbers. Calc has no contour chart; the coloured grid shows the same information.
- **3D surface:** the same grid drawn as a surface (section 6.3), inserted as a picture.

Factors not on the map are held at their centre (coded 0); the map title says at which values.

#### Checked against

Montgomery's Example 11.2 (CCD, time and temperature, yield):
- fitted model 79.94 + 0.995x₁ + 0.515x₂ − 1.376x₁² − 1.001x₂² + 0.250x₁x₂;
- lack of fit SS 0.2833 (3 df), pure error SS 0.2120 (4 df);
- stationary point at 86.95 min and 176.53 °F with predicted yield 80.21, eigenvalues −0.9634 and −1.4141 (a maximum);
- equation in actual units −1430.52 + 7.81ξ₁ + 13.27ξ₂ − 0.055ξ₁² − 0.040ξ₂² + 0.010ξ₁ξ₂.

All of these are reproduced from a Design sheet made by LibreStats.

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
| `qq` | `GroupRange` (empty), `Positions` = `blom`/`vdw`/`tukey`/`rankit` (blom), `Line` = `quartiles`/`meansd`/`none` (quartiles), `Title`, `XTitle`, `YTitle` (empty = automatic), `Color` (RGB number) |
| `surface` | Range: three columns x, y, z. `Azimuth` (−60), `Elevation` (30), `MeshLines` (True), `Title`, `XTitle`, `YTitle`, `ZTitle` (empty = automatic) |
| `normality` | `GroupRange` (empty), `Alpha` (0.05) |
| `variances` | `GroupRange` (empty), `Center` = `median`/`mean` (median), `Alpha` (0.05) |
| `ttestOne` | `GroupRange` (empty), `Mu` (0), `Alternative` = `two-sided`/`greater`/`less` (two-sided), `Alpha` (0.05) |
| `ttestWelch` | `GroupRange`, `Alternative`, `Alpha` |
| `ttestPaired` | `Alternative`, `Alpha` |
| `anovaOneWay`, `anovaWelch` | `GroupRange`, `Alpha` |
| `anovaRM` | `Alpha` |
| `posthocTukey`, `posthocGamesHowell` | `GroupRange`, `Alpha` |
| `posthocDunnett` | `GroupRange`, `Control` (empty = first group), `Alternative`, `Alpha` |
| `posthocDunn` | `GroupRange`, `Adjust` = `holm`/`bonferroni`/`sidak`/`bh` (holm), `Alpha` |
| `doeDesign` | Range: the factor table. `Design` = `full`/`fractional`/`pb`/`ccd`/`bbd` (full), `Runs` (empty = automatic), `CentrePoints` (screening: 0; response surface: empty = automatic), `Replicates` (1), `Blocks` (1), `Randomize` (True), `Seed` (empty = random), `Responses` (comma-separated, "Response"). CCD: `Alpha` = `rotatable`/`orthogonal`/`face`/a number (rotatable), `FractionalCube` (True), `Inscribed` (False), `AxialCentrePoints` (empty = automatic) |
| `doeAnalyseScreening` | Range: the design table (header row included). `Response` (empty = first response column), `Model` = `2fi`/`main` (2fi), `Alpha` (0.05), `Charts` (True) |
| `doeAnalyseRSM` | Range: the design table (header row included). `Response` (empty = first response column), `Model` = `quadratic`/`interaction`/`linear` (quadratic), `Alpha` (0.05), `XFactor`, `YFactor` (factor names or letters for the response map; empty = first two), `Charts` (True) |

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

## 11. Limitations

- **Normality assumptions:** Grubbs, Dixon, Mahalanobis (χ² cutoff), Little's test and MICE's regression models assume (multivariate) normality; this is stated in the reports. Consider transforming skewed data first (section 4.1).
- **Mahalanobis distance** uses the classical mean and covariance, which outliers themselves can distort (masking). A robust version (MCD; Leys et al., 2018) is not yet available.
- **No pooling of MICE results yet:** LibreStats creates the imputed datasets but cannot yet combine analyses with Rubin's rules.
- **MAR vs MNAR** cannot be distinguished by Little's test or any other test on the observed data.
- **Performance:** pure Python is slower than compiled libraries. Data-preparation tools handle typical datasets (thousands of rows) easily; MICE with many variables, imputations and iterations takes longer.
- **DOE / RSM is new and still being tested.** The designs and analyses are checked against published examples (NIST, Montgomery, Box & Hunter, Box & Behnken) and statsmodels, but the module has had little use on real experiments so far. Check the design summary before running an experiment, and report anything unexpected.
- **DOE / RSM scope:** two-level screening designs and central composite / Box-Behnken designs only (no definitive screening, Doehlert, mixture or optimal designs yet). Factors with text levels are possible in screening designs only. Each response is analysed on its own: there is no simultaneous optimisation of several responses (desirability) yet. The response map and 3D surface show two factors at a time, with the others at their centre.
- **Platforms:** developed and tested on Linux with LibreOffice 26.8; not yet tested on Windows or macOS.

---

## References

- Aitchison, J. (1982). The statistical analysis of compositional data. *Journal of the Royal Statistical Society B*, 44(2), 139–177.
- Bakeman, R. (2005). Recommended effect size statistics for repeated measures designs. *Behavior Research Methods*, 37(3), 379–384.
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B*, 57(1), 289–300.
- Blom, G. (1958). *Statistical Estimates and Transformed Beta-Variables*. Wiley.
- Bowyer, A. (1981). Computing Dirichlet tessellations. *The Computer Journal*, 24(2), 162–166.
- Box, G. E. P., & Behnken, D. W. (1960). Some new three level designs for the study of quantitative variables. *Technometrics*, 2(4), 455–475.
- Box, G. E. P., & Cox, D. R. (1964). An analysis of transformations. *Journal of the Royal Statistical Society B*, 26(2), 211–252.
- Box, G. E. P., & Hunter, J. S. (1957). Multi-factor experimental designs for exploring response surfaces. *Annals of Mathematical Statistics*, 28(1), 195–241.
- Box, G. E. P., & Hunter, J. S. (1961). The 2^(k−p) fractional factorial designs. *Technometrics*, 3(3), 311–351 and 3(4), 449–458.
- Box, G. E. P., & Wilson, K. B. (1951). On the experimental attainment of optimum conditions. *Journal of the Royal Statistical Society B*, 13(1), 1–45.
- Brown, M. B., & Forsythe, A. B. (1974). Robust tests for the equality of variances. *Journal of the American Statistical Association*, 69(346), 364–367.
- Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences* (2nd ed.). Erlbaum.
- Daniel, C. (1959). Use of half-normal plots in interpreting factorial two-level experiments. *Technometrics*, 1(4), 311–341.
- Delacre, M., Lakens, D., & Leys, C. (2017). Why psychologists should by default use Welch's t-test instead of Student's t-test. *International Review of Social Psychology*, 30(1), 92–101.
- Delacre, M., Leys, C., Mora, Y. L., & Lakens, D. (2019). Taking parametric assumptions seriously: Arguments for the use of Welch's F-test instead of the classical F-test in one-way ANOVA. *International Review of Social Psychology*, 32(1), 13.
- Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society B*, 39(1), 1–38.
- Dixon, W. J. (1951). Ratios involving extreme values. *Annals of Mathematical Statistics*, 22(1), 68–78.
- Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics*, 6(3), 241–252.
- Dunnett, C. W. (1955). A multiple comparison procedure for comparing several treatments with a control. *Journal of the American Statistical Association*, 50(272), 1096–1121.
- Dunnett, C. W. (1964). New tables for multiple comparisons with a control. *Biometrics*, 20(3), 482–491.
- Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal, C. (2003). Isometric logratio transformations for compositional data analysis. *Mathematical Geology*, 35(3), 279–300.
- Filzmoser, P., Hron, K., & Templ, M. (2018). *Applied Compositional Data Analysis*. Springer.
- Fisher, R. A. (1925). *Statistical Methods for Research Workers*. Oliver & Boyd.
- Fries, A., & Hunter, W. G. (1980). Minimum aberration 2^(k−p) designs. *Technometrics*, 22(4), 601–608.
- Games, P. A., & Howell, J. F. (1976). Pairwise multiple comparison procedures with unequal n's and/or variances: A Monte Carlo study. *Journal of Educational Statistics*, 1(2), 113–125.
- Girden, E. R. (1992). *ANOVA: Repeated Measures*. Sage.
- Greenhouse, S. W., & Geisser, S. (1959). On methods in the analysis of profile data. *Psychometrika*, 24(2), 95–112.
- Grubbs, F. E. (1969). Procedures for detecting outlying observations in samples. *Technometrics*, 11(1), 1–21.
- Hoaglin, D. C., & Iglewicz, B. (1987). Fine-tuning some resistant rules for outlier labeling. *Journal of the American Statistical Association*, 82(400), 1147–1149.
- Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics*, 6(2), 65–70.
- Huynh, H., & Feldt, L. S. (1976). Estimation of the Box correction for degrees of freedom from sample data in randomized block and split-plot designs. *Journal of Educational Statistics*, 1(1), 69–82.
- Kramer, C. Y. (1956). Extension of multiple range tests to group means with unequal numbers of replications. *Biometrics*, 12(3), 307–310.
- Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association*, 47(260), 583–621.
- Lakens, D. (2013). Calculating and reporting effect sizes to facilitate cumulative science: A practical primer for t-tests and ANOVAs. *Frontiers in Psychology*, 4, 863.
- Lecoutre, B. (1991). A correction for the ε̃ approximate test in repeated measures designs with two or more independent groups. *Journal of Educational Statistics*, 16(4), 371–372.
- Lenth, R. V. (1989). Quick and easy analysis of unreplicated factorials. *Technometrics*, 31(4), 469–473.
- Levene, H. (1960). Robust tests for equality of variances. In I. Olkin (Ed.), *Contributions to Probability and Statistics* (pp. 278–292). Stanford University Press.
- Leys, C., Klein, O., Dominicy, Y., & Ley, C. (2018). Detecting multivariate outliers: Use a robust variant of the Mahalanobis distance. *Journal of Experimental Social Psychology*, 74, 150–156.
- Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. *Journal of Experimental Social Psychology*, 49(4), 764–766.
- Little, R. J. A. (1988). A test of missing completely at random for multivariate data with missing values. *Journal of the American Statistical Association*, 83(404), 1198–1202.
- Mahalanobis, P. C. (1936). On the generalised distance in statistics. *Proceedings of the National Institute of Sciences of India*, 2, 49–55.
- Martín-Fernández, J. A., Barceló-Vidal, C., & Pawlowsky-Glahn, V. (2003). Dealing with zeros and missing values in compositional data sets using nonparametric imputation. *Mathematical Geology*, 35(3), 253–278.
- Mauchly, J. W. (1940). Significance test for sphericity of a normal n-variate distribution. *Annals of Mathematical Statistics*, 11(2), 204–209.
- Montgomery, D. C. (2017). *Design and Analysis of Experiments* (9th ed.). Wiley.
- Myers, R. H., Montgomery, D. C., & Anderson-Cook, C. M. (2016). *Response Surface Methodology: Process and Product Optimization Using Designed Experiments* (4th ed.). Wiley.
- Newell, M. E., Newell, R. G., & Sancha, T. L. (1972). A solution to the hidden surface problem. *Proceedings of the ACM Annual Conference*, 1, 443–450.
- NIST/SEMATECH (2012). *e-Handbook of Statistical Methods*, section 5.3.3.4.7: Summary tables of useful fractional factorial designs. https://www.itl.nist.gov/div898/handbook/
- Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared statistics: Measures of effect size for some common research designs. *Psychological Methods*, 8(4), 434–447.
- Plackett, R. L., & Burman, J. P. (1946). The design of optimum multifactorial experiments. *Biometrika*, 33(4), 305–325.
- Royston, P. (1992). Approximating the Shapiro-Wilk W-test for non-normality. *Statistics and Computing*, 2(3), 117–119.
- Royston, P. (1995). Remark AS R94: A remark on algorithm AS 181: The W-test for normality. *Applied Statistics*, 44(4), 547–551.
- Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.
- Satterthwaite, F. E. (1946). An approximate distribution of estimates of variance components. *Biometrics Bulletin*, 2(6), 110–114.
- Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). *Biometrika*, 52(3–4), 591–611.
- Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. *Journal of the American Statistical Association*, 62(318), 626–633.
- Smithson, M., & Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. *Psychological Methods*, 11(1), 54–71.
- Student. (1908). The probable error of a mean. *Biometrika*, 6(1), 1–25.
- Sturges, H. A. (1926). The choice of a class interval. *Journal of the American Statistical Association*, 21(153), 65–66.
- Tabachnick, B. G., & Fidell, L. S. (2013). *Using Multivariate Statistics* (6th ed.). Pearson.
- Troyanskaya, O., Cantor, M., Sherlock, G., Brown, P., Hastie, T., Tibshirani, R., Botstein, D., & Altman, R. B. (2001). Missing value estimation methods for DNA microarrays. *Bioinformatics*, 17(6), 520–525.
- Tukey, J. W. (1953). *The Problem of Multiple Comparisons*. Unpublished manuscript, Princeton University.
- Tukey, J. W. (1977). *Exploratory Data Analysis*. Addison-Wesley.
- van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate imputation by chained equations in R. *Journal of Statistical Software*, 45(3), 1–67.
- van den Berg, R. A., Hoefsloot, H. C. J., Westerhuis, J. A., Smilde, A. K., & van der Werf, M. J. (2006). Centering, scaling, and transformations: Improving the biological information content of metabolomics data. *BMC Genomics*, 7, 142.
- van der Walt, S., & Smith, N. (2015). *mpl colormaps* (viridis). https://bids.github.io/colormap/
- Wang, H., & Song, M. (2011). Ckmeans.1d.dp: Optimal k-means clustering in one dimension by dynamic programming. *The R Journal*, 3(2), 29–33.
- Watson, D. F. (1981). Computing the n-dimensional Delaunay tessellation with application to Voronoi polytopes. *The Computer Journal*, 24(2), 167–172.
- Welch, B. L. (1947). The generalization of "Student's" problem when several different population variances are involved. *Biometrika*, 34(1–2), 28–35.
- Welch, B. L. (1951). On the comparison of several mean values: An alternative approach. *Biometrika*, 38(3–4), 330–336.
- Yeo, I.-K., & Johnson, R. A. (2000). A new family of power transformations to improve normality or symmetry. *Biometrika*, 87(4), 954–959.
