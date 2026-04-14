# ROCsurvcomp

## Overview

**ROCsurvcomp** is an R package for ROC-based methods to compare survival curves under non-proportional hazards (non-PH).
Traditional log-rank test may lose power when the proportional hazards assumption is violated. Other methods, such as the Flemming-Harrington (FH) family of weighted log-rank tests require prior knowledge of the underlying non-PH patterns, and mis-specified patterns may lead to a substantial loss of statistical power. This package provides alternative approaches for comparing survival curves based on ROC curve, without requiring for any prior knowledge of the underlying non-PH pattern.

## Key Features

- Implements nonparametric kernel-based approaches for comparing two survival distributions
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
# Example (replace with your actual data)

result <- surv.comp(
  time = survival_time,
  status = censoring_status,
  group = group_indicator,
  n_perm = 1000,
  censor_type = "double",
  method = "joint_method")

print(result)
```

## Methodology

This package implements ROC-based approaches for comparing survival curves:

- **ROC Length**: Measures separation between survival distributions
- **Overlap Coefficient (OVL)**: Quantifies similarity between distributions
- **Joint Test**: Combines ROC length and OVL for improved inference

These methods are especially useful when:
- Proportional hazards assumption is violated
- Treatment effects are not constant over time
- Survival curves cross


---

## Authors

**Mohammod Mahmudur Rahman**<br>
PhD Student, Department of Biostatistics & Data Science<br>
University of Kansas Medical Center

**Leonidas Bantis**<br>
Associate Professor, Department of Biostatistics & Data Science<br>
University of Kansas Medical Center


## Status

🚧 This package is under active development.

