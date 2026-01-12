# Histogram-Based Mutual Information (MATLAB & R)

This repository contains **matched MATLAB and R implementations** of a histogram-based estimator of **mutual information (MI)** for continuous variables.

The goal is to provide a **transparent, reproducible, and cross-language-consistent** implementation that:

* uses identical binning rules,
* applies Laplace smoothing (Jeffrey-Perks law),
* and yields comparable (albein not perfectly identical) numerical results in MATLAB and R.

---


## Overview

Mutual information measures statistical dependence between two variables:

<img width="451" height="94" alt="Screenshot 2026-01-12 at 11 40 13" src="https://github.com/user-attachments/assets/4923a900-9aee-4e54-ac08-d00c3fa261aa" />

In this repository:

* probabilities are estimated via a **2D histogram**,
* **Laplace smoothing is applied to the joint histogram only**,
* marginals are derived from the smoothed joint distribution.

---

## Key Implementation Choices

### 1. Histogram binning

* Bin edges are **equally spaced between min and max**, using:

  * MATLAB: `linspace(min, max, nBins+1)`
  * R: `seq(min, max, length.out = nBins+1)`
* Bin intervals follow MATLAB’s convention:

  * **left-closed, right-open** `[e_i, e_{i+1})`
  * final bin includes the right edge

### 2. Smoothing

* **Laplace smoothing** is applied to the **joint histogram only**:

  ```text
  pXY = (counts + smoothingValue) / sum(counts + smoothingValue)
  ```
* Marginals are computed from `pXY`, not smoothed independently.

### 3. Orientation consistency

* Joint histogram is treated as **(Y bins × X bins)** in both languages.
* This matches MATLAB’s `histogram2.Values` layout.

### 4. Units

* `"bits"` → log base 2
* `"nats"` → natural logarithm

---

## Files

```
.
├── getMI.m        # MATLAB implementation
├── getMI.R        # R implementation
├── x1.csv         # Example independent data
├── y1.csv
├── x2.csv         # Example dependent data
├── y2.csv
└── README.md
```

---

## MATLAB Usage

### Function

```matlab
MI = getMI(x, y, nBins, smoothingValue, units)
```

### Example

```matlab
x = randn(1000,1);
y = randn(1000,1);
MI = getMI(x, y, 20, 0.5, "bits");
```

### Monte-Carlo sanity check

```matlab
nBins = 25;
smoothingValue = 0.5;
nRep = 100;
nObs = 10000;

for k = 1:nRep
    x = randn(nObs,1);
    y = randn(nObs,1);
    MI(k) = getMI(x,y,nBins,smoothingValue);
end
histogram(MI)
```

---

## R Usage

### Function

```r
MI <- getMI(x, y, nBins = 20, smoothingValue = 0.5, units = "bits")
```

### Example

```r
set.seed(1310)
x <- rnorm(1000)
y <- rnorm(1000)

source("getMI.R")
MI_bits <- getMI(x, y, nBins = 20, smoothingValue = 0.5, units = "bits")
```

### Using pre-generated MATLAB data

```r
x1 <- as.vector(read.csv("x1.csv", header = FALSE))
y1 <- as.vector(read.csv("y1.csv", header = FALSE))
x2 <- as.vector(read.csv("x2.csv", header = FALSE))
y2 <- as.vector(read.csv("y2.csv", header = FALSE))

MI1 <- getMI(x1, y1, nBins = 20, smoothingValue = 0.5, units = "bits")
MI2 <- getMI(x2, y2, nBins = 20, smoothingValue = 0.5, units = "bits")
```

---

## Important Notes & Limitations

### Histogram MI is **bin-sensitive**

* Results depend strongly on:

  * number of bins,
  * bin placement,
  * data range.
* Heavy-tailed variables (e.g. `y = x + 1/x`) can collapse into few bins when using min/max binning, yielding **low MI even for dependent variables**.

### Finite-sample bias

* Histogram MI is **positively biased** for small sample sizes.
* Independent variables typically yield small but nonzero MI.

### This is not a kNN estimator

If you need:

* bin-free estimation,
* better behavior for nonlinear or heavy-tailed dependencies,

consider **KSG / k-nearest-neighbor MI estimators** instead.

---

## Purpose

This code is intended for:

* didactic use,
* controlled simulations,
* cross-language reproducibility (MATLAB ↔ R),
* understanding estimator behavior and pitfalls.

It is **not** meant to be a state-of-the-art MI estimator.

---

## References

Burroni, F., & Tilsen, S. (2025). Thai speakers time lexical tones to supralaryngeal articulatory events. Journal of Phonetics, 108, 101389. https://doi.org/https://doi.org/10.1016/j.wocn.2024.101389 

Iskarous, K., Mooshammer, C., Hoole, P., Recasens, D., Shadle, C. H., Saltzman, E., & Whalen, D. H. (2013). The coarticulation/invariance scale: Mutual information as a measure of coarticulation resistance, motor synergy, and articulatory invariance. The Journal of the Acoustical Society of America, 134(2), 1271-1282. 


---

## License

MIT

