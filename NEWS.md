# ROCsurvcomp 0.1.2
## Changes suggested by CRAN
- Updated the Description field to expand the acronym ROC to Receiver Operating Characteristic (ROC).
- Removed `set.seed()` from three functions (`doublecen_survivaltest_functions.R`, `leftcen_survivaltest_functions.R`, `rightcen_survivaltest_functions.R`).

## Additional Minor changes
- Added additional references to the vignette.
- Added an explanatory text to the examples of the `surv.comp()` function in the help file.


# ROCsurvcomp 0.1.1
## Minor changes
- Updated license field to use GPL-3 without an additional LICENSE file
- Added references with DOIs in the Description field

# ROCsurvcomp 0.1.0
## Initial release
- Provides ROC-based methods for comparing survival distributions under non-proportional hazards
- Implements ROC length, overlap coefficient (OVL), and joint ROC length–OVL tests
- Supports right, left, and doubly censored data
- Includes permutation-based inference procedures
- Contains vignette with examples and methodological details
