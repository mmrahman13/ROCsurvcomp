######################################################################
####### Kernel Density Estimation with Right censoring ########
######################################################################

kde_right_censored <- function(x, data, censor) {
  # x: Points at which to estimate the density
  # data: Observed time-to-event data
  # censoring: 0 if event, 1 if censored

  # Step 1: Compute the Kaplan-Meier estimator
  fit <- survival::survfit(survival::Surv(data, censor == 0) ~ 1)
  F_km <- 1 - fit$surv                         # Convert survival to CDF
  F_km_jumps <- c(F_km[1], diff(F_km))         # KM jumps at each event time
  data_uncen <- fit$time                       # Event times where KM jumps

  # Step 2: Bandwidth selection using Silverman's rule of thumb for bandwidth
  n <- length(data)                            # Number of TOTAL observations
  k <- 1 / qnorm(3/4)
  median_index <- which(F_km >= 0.50)[1]       # Find first time where F_km >= 0.50
  median_data_uncen <- data_uncen[median_index]

  sigma_hat <- k * median(abs(data - median_data_uncen))    # Median Absolute Deviation (MAD)
  h <- (4 / (3 * n))^(1/5) * sigma_hat

  # Step 3: Compute the weighted kernel density estimate
  pdf <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    kernel_values <- (1 / (sqrt(2 * pi) * h)) * exp(-0.5 * ((xi - data_uncen) / h)^2)
    pdf[i] <- sum(F_km_jumps * kernel_values)               # Weighted sum using KM jumps
  }

  return(pdf)
}




######################################################################
####### Kernel Density Estimation with Double censoring ########
######################################################################

kde_double_censored <- function(x, data, censor) {
  # x: Points at which to estimate the density
  # data: Observed time-to-event data
  # censoring: 0 if event, 1 if right censored, -1 if left censored

  # Initialize offset to 0
  offset <- 0

  # Step 0: Shift data ONLY if it has negative obs (to satisfy "icfit" function)
  min_obs <- min(data)
  if (min_obs < 0) {
    offset <- abs(min_obs) + 1
  }

  data_shifted <- data + offset

  # Step 1: Compute the Turnbull estimator using shifted data
  lower <- ifelse(censor == -1, -Inf, data_shifted)
  upper <- ifelse(censor == -1, data_shifted,
                  ifelse(censor == 0, data_shifted,
                         ifelse(censor == 1, Inf, NA)))

  fit <- interval::icfit(survival::Surv(lower, upper, type = "interval2") ~ 1)
  F_tbull <- as.numeric(cumsum(fit$pf))                 # CDF
  F_tbull_jumps <- c(F_tbull[1], diff(F_tbull))         # Turnbull jumps at each event time

  # Extract jump times and shift them back to original scale
  data_uncen_shifted <- as.numeric(sub(".*,", "", gsub("\\[|\\]", "", names(cumsum(fit$pf)))))
  data_uncen <- data_uncen_shifted - offset

  # Step 2: Bandwidth selection using Silverman's rule of thumb for bandwidth
  n <- length(data)
  k <- 1 / qnorm(3/4)
  median_index <- which(F_tbull >= 0.50)[1]       # Find first time where F_tbull >= 0.50
  median_data_uncen <- data_uncen[median_index]

  sigma_hat <- k * median(abs(data - median_data_uncen))
  h <- (4 / (3 * n))^(1/5) * sigma_hat

  # Step 3: Compute the weighted kernel density estimate
  pdf <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i] # This is the original x provided by user
    kernel_values <- (1 / (h * sqrt(2 * pi))) * exp(-0.5 * ((xi - data_uncen) / h)^2)
    pdf[i] <- sum(F_tbull_jumps * kernel_values)
  }

  return(pdf)
}


