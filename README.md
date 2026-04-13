# ROCsurvcomp

## Overview

**ROCsurvcomp** is an R package for ROC-based methods to compare survival curves under non-proportional hazards (non-PH).
Traditional methods such as the log-rank test may lose power when the proportional hazards assumption is violated. This package provides alternative approaches based on ROC curve characteristics.

## Key Features

- ROC-based comparison of survival curves
- Designed for *non-proportional hazards* settings
- Implements:
  - ROC length-based test
  - Overlap (OVL) measure
  - Joint ROC–OVL testing framework
- Supports right, left and doubly-censored survival data
- Permutation-based inference


## Installation

You can install from GitHub:

```r
# install.packages("devtools")    # if not installed
devtools::install_github("mmrahman13/ROCsurvcomp")
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
- Survival curves cross
- Treatment effects vary over time
- Proportional hazards assumption is violated

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


## Citation

If you use this package in your research, please cite appropriately (details to be added).
