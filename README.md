# ROCsurvcomp

## Overview

**ROCsurvcomp** is an R package for comparing survival curves under non-proportional hazards (non-PH) using ROC-based methods.
Traditional log-rank test may lose power when the proportional hazards assumption is violated. Other methods, such as the Flemming-Harrington (FH) family of weighted log-rank tests require prior knowledge of the underlying non-PH patterns, and mis-specified patterns may lead to a substantial loss of statistical power. This package provides alternative approaches for comparing survival curves based on ROC curve, without requiring prior knowledge of the underlying non-PH pattern, and can accommodate right, left, and doubly censored data.

## Key Features

- Implements nonparametric and semiparametric kernel-based approaches for comparing two survival distributions
- Designed for handling *non-proportional hazards* settings
- Methods include:
  - ROC length-based test
  - Overlap coefficient (OVL)-based test
  - Joint ROC length-OVL-based test
- Supports right, left and doubly-censored survival data
- Permutation-based inference with two-sided hypothesis testing


## Installation

You can install from GitHub:

```r
# install.packages("remotes")    # if not installed
remotes::install_github("mmrahman13/ROCsurvcomp")
```


## Usage

### Basic Example

```r
# Example

library(ROCsurvcomp)
data("EarlyEffectData")

surv.comp(time = EarlyEffectData$time,
          status = EarlyEffectData$status,
          group = EarlyEffectData$group,
          censor_type = "right",
          method = "roc_length",
          n_perm = 1000,
          progress = TRUE,
          plot = FALSE)
```

## Methodology

This package implements ROC-based approaches for comparing two survival curves:

- **ROC Length**: Measures separation between two distributions without relying on any stochastic ordering
- **Overlap Coefficient (OVL)**: Quantifies the similarity between two distributions using the common area between two probability density functions
- **Joint Test**: Combines ROC length and OVL-based methods by constructing convex hull of the permuted samples based on their Euclidean distance from the origin.

These methods are especially useful when:
- The proportional hazards assumption is violated
- Treatment effects are not constant over time and no prior information about the pattern of effects is available
- Survival curves cross.


---

## Authors

**Mohammod Mahmudur Rahman**<br>
PhD Student, Department of Biostatistics & Data Science<br>
University of Kansas Medical Center

**Leonidas Bantis**<br>
Associate Professor, Department of Biostatistics & Data Science<br>
University of Kansas Medical Center

