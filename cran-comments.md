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
- Changed the License to GPL-3 and removed the LICENSE file.
- Added references with DOIs for the Bantis (2021) and Franco-Pereira (2021) papers to the Description field.

