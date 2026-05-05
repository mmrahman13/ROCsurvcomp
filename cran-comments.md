## Test environments
- macOS 15.7.2 (Apple Silicon), R 4.5.1
- Local machine

## R CMD check results
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Package description
The ROCsurvcomp package provides nonparametric and semiparametric methods for comparing survival distributions under non-proportional hazards. It implements ROC length, overlap coefficient (OVL), and a joint ROC length–OVL-based inference framework, supporting right, left, and doubly censored data.

## Additional comments
- The vignette builds successfully and includes all required dependencies.
- A small number of permutations (n_perm = 10) is used in all examples and vignettes to ensure that all code runs within CRAN's time limits without requiring \donttest{} tags.
- No external data, APIs, or system-specific features are used.

## Resubmission
This is a resubmission. Following the maintainer's suggestions, I have:
- Updated the Description field to expand the acronym ROC to Receiver Operating Characteristic (ROC).
- Removed `set.seed()` from three functions (`doublecen_survivaltest_functions.R`, `leftcen_survivaltest_functions.R`, `rightcen_survivaltest_functions.R`).

I have also performed some additional minor changes:
- Added additional references to the vignette.
- Added an explanatory text to the examples of the `surv.comp()` function in the help file.
