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
library(PWEXP)

# Generating data with crossing survivals
set.seed(800)
n_trt <- 100
break_trt <- c(3, 7)
rate_trt <- c(log(2)/5, log(2)/8, log(2)/13)
rate.censor_trt <- c(log(2)/20.80, log(2)/23.80, log(2)/28.80)
event_trt <- PWEXP::rpwexp(n_trt, rate = rate_trt, breakpoint = break_trt)
censor_trt <- PWEXP::rpwexp(n_trt, rate = rate.censor_trt, breakpoint = break_trt)

n_ctrl <- 100
rate_ctrl <- log(2)/8
rate.censor_ctrl <- log(2)/24
event_ctrl <- rexp(n_ctrl, rate = rate_ctrl)
censor_ctrl <- rexp(n_ctrl, rate = rate.censor_ctrl)

# Observed time and censoring status (0 = event, 1 = right-censored)
time_trt <- pmin(event_trt, censor_trt)
status_trt <- ifelse(event_trt <= censor_trt, 0, 1)
time_ctrl <- pmin(event_ctrl, censor_ctrl)
status_ctrl <- ifelse(event_ctrl <= censor_ctrl, 0, 1)
time <- c(time_trt, time_ctrl)
status <- c(status_trt, status_ctrl)
group <- c(rep(1, n_trt), rep(2, n_ctrl))

# Run `surv.comp` function
surv.comp(
  time = time,
  status = status,
  group = group,
  censor_type = "right",
  method = "joint_method",
  n_perm = 10,
  progress = TRUE,
  plot = TRUE
)
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

