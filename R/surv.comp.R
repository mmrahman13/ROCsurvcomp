#' @title ROC-Based Methods for Comparing Survival Distributions with Right, Left, and Doubly Censored Data
#'
#' @description Performs a nonparametric comparison of two survival distributions under
#' right, left, or double censoring using ROC-based metrics, including ROC curve length,
#' overlap coefficient (OVL), and a joint ROC length–OVL test.
#'
#' @importFrom stats dnorm integrate median optim pnorm qnorm rnorm sd
#' @importFrom utils tail
#'
#' @param time Numeric vector of observed follow-up times (event or censoring times)
#'   for all observations.
#' @param status Numeric vector indicating censoring status for each observation:
#'   \itemize{
#'     \item For \code{censor_type = "right"}: 0 = event, 1 = right-censored
#'     \item For \code{censor_type = "left"}: 0 = event, -1 = left-censored
#'     \item For \code{censor_type = "double"}: 0 = event, 1 = right-censored, -1 = left-censored.
#'   }
#' @param group Numeric vector indicating the group label for each observation.
#'   Must contain exactly two groups coded as 1 and 2.
#' @param censor_type Character string specifying the type of censoring.
#'   Must be either \code{"right"}, \code{"left"}, or \code{"double"}.
#' @param method Character string specifying the test to perform.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"roc_length"}: ROC curve length-based test
#'     \item \code{"ovl"}: overlap coefficient-based test
#'     \item \code{"joint_method"}: joint ROC length and OVL test.
#'   }
#' @param n_perm Integer specifying the number of permutation samples used to compute p-values.
#' @param progress Logical value indicating whether to display a progress bar
#'   during the permutation test. Default is \code{TRUE}. If \code{FALSE},
#'   the computation runs silently without showing progress updates.
#' @param plot Logical value indicating whether to generate a convex hull plot
#'   based on the joint ROC length–OVL test. Default is \code{FALSE}.
#'   Plotting is only available when \code{method = "joint_method"}; otherwise,
#'   an error is returned if \code{plot = TRUE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{message}: A character string describing the testing procedure.
#'   \item \code{result}: A data frame with rows corresponding to the methods
#'   (ROC Length, OVL, and/or Joint ROC Length-OVL) and columns:
#'   \itemize{
#'     \item \code{estimate}: Estimated values of ROC length and/or OVL
#'     \item \code{p_value}: Permutation-based two-sided p-values.
#'   }
#'   \item \code{plot}: A \code{ggplot2} object showing the convex hull
#'   visualization of the permutation distribution for the joint ROC length-OVL test.
#'   This is returned only when \code{method = "joint_method"} and
#'   \code{plot = TRUE}; otherwise, it is \code{NULL}.
#' }
#'
#' @details
#' This function implements permutation-based two-sided hypothesis tests for comparing
#' two survival distributions without relying on proportional hazards assumptions.
#' The ROC length measures global separation between survival curves, while the
#' overlap coefficient (OVL) quantifies distributional similarity. Density estimation is performed using nonparametric kernel methods,
#' where the cumulative distribution function (CDF) is estimated
#' using the Kaplan–Meier estimator under right censoring or the Turnbull estimator under left and double censoring.
#'
#' Input validation is performed to ensure consistency of vector lengths,
#' censoring indicators, and method specification.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data("EarlyEffectData")
#' surv.comp(time = EarlyEffectData$time, status = EarlyEffectData$status,
#'           group = EarlyEffectData$group, censor_type = "right",
#'           method = "roc_length", n_perm = 1000, progress = TRUE, plot = FALSE)
#' }
#' @export
surv.comp <- function(time, status, group, censor_type, method,
                      n_perm, progress = TRUE, plot = FALSE) {

  ## -----------------------------
  ## Check required inputs
  ## -----------------------------
  if (missing(time) || missing(status) || missing(group) || missing(n_perm) ||
      missing(censor_type) || missing(method)) {
    stop("All arguments ('time', 'status', 'group', 'n_perm', 'censor_type', and 'method') must be provided.")
  }


  if (censor_type == "right") {

    ## -----------------------------
    ## 1. Check lengths
    ## -----------------------------
    if (length(unique(c(length(time), length(status), length(group)))) != 1) {
      stop("Length mismatch: 'time', 'status', and 'group' must have the same length.")
    }

    ## -----------------------------
    ## 2. Check missing values
    ## -----------------------------
    if (any(is.na(c(time, status, group)))) {
      stop("Missing values found in at least one of 'time', 'status', or 'group'.")
    }

    ## -----------------------------
    ## 3. Check data types
    ## -----------------------------
    if (!is.numeric(time)) {
      stop("'time' must be a numeric vector.")
    }

    if (!is.numeric(status)) {
      stop("'status' must be a numeric vector.")
    }

    if (!is.numeric(group)) {
      stop("'group' must be a numeric vector.")
    }


    ## -----------------------------
    ## 4. Check censoring values
    ## -----------------------------
    if (!all(status %in% c(0, 1))) {
      stop("'status' must contain only 0 (event) and 1 (right-censored) values, when censor_type = 'right'.")
    }

    ## -----------------------------
    ## 5. Check n_perm
    ## -----------------------------
    if (!is.numeric(n_perm) || length(n_perm) != 1 || n_perm <= 0 || n_perm %% 1 != 0) {
      stop("'n_perm' must be a positive integer.")
    }

    ## -----------------------------
    ## 6. Check method
    ## -----------------------------
    valid_methods <- c("roc_length", "ovl", "joint_method")

    if (length(method) != 1 || !method %in% valid_methods) {
      stop("'method' must be exactly one of: 'roc_length', 'ovl', or 'joint_method'")
    }

    ## -----------------------------
    ## 7. Check "group" values
    ##    & filter data based on that
    ## -----------------------------
    if (!all(group %in% c(1, 2)) || !all(c(1, 2) %in% group)) {
      stop("'group' must contain exactly two groups coded as 1 and 2.")
    }

    ## -----------------------------
    ## 8. Check progress bar
    ## -----------------------------
    if (!is.logical(progress) || length(progress) != 1) {
      stop("'progress' must be either TRUE or FALSE.")
    }

    ## -----------------------------
    ## 9. Check plot argument
    ## -----------------------------
    if (!is.logical(plot) || length(plot) != 1) {
      stop("'plot' must be either TRUE or FALSE.")
    }

    if (plot && method != "joint_method") {
      stop("Plotting is only available for method = 'joint_method'. Set plot = FALSE or change method.")
    }

    time1 <- time[group == 1]
    censor1 <- status[group == 1]
    time2 <- time[group == 2]
    censor2 <- status[group == 2]

    ## -----------------------------
    ## Run selected method
    ## -----------------------------
    if (method == "roc_length") {
      return(roc_RightCenSurvival_test(time1, censor1, time2, censor2, n_perm, progress))
    }

    if (method == "ovl") {
      return(ovl_RightCenSurvival_test(time1, censor1, time2, censor2, n_perm, progress))
    }

    if (method == "joint_method") {
      return(joint.roc_ovl_RightCenSurvival_test(time1, censor1, time2, censor2, n_perm, progress, plot))
    }
  }


  else if (censor_type == "left") {

    ## -----------------------------
    ## 1. Check lengths
    ## -----------------------------
    if (length(unique(c(length(time), length(status), length(group)))) != 1) {
      stop("Length mismatch: 'time', 'status', and 'group' must have the same length.")
    }

    ## -----------------------------
    ## 2. Check missing values
    ## -----------------------------
    if (any(is.na(c(time, status, group)))) {
      stop("Missing values found in at least one of 'time', 'status', or 'group'.")
    }

    ## -----------------------------
    ## 3. Check data types
    ## -----------------------------
    if (!is.numeric(time)) {
      stop("'time' must be a numeric vector.")
    }

    if (!is.numeric(status)) {
      stop("'status' must be a numeric vector.")
    }

    if (!is.numeric(group)) {
      stop("'group' must be a numeric vector.")
    }

    ## -----------------------------
    ## 4. Check censoring values
    ## -----------------------------
    if (!all(status %in% c(0, -1))) {
      stop("'status' must contain only 0 (event) and -1 (left-censored) values, when censor_type = 'left'.")
    }

    ## -----------------------------
    ## 5. Check n_perm
    ## -----------------------------
    if (!is.numeric(n_perm) || length(n_perm) != 1 || n_perm <= 0 || n_perm %% 1 != 0) {
      stop("'n_perm' must be a positive integer.")
    }

    ## -----------------------------
    ## 6. Check method
    ## -----------------------------
    valid_methods <- c("roc_length", "ovl", "joint_method")

    if (length(method) != 1 || !method %in% valid_methods) {
      stop("'method' must be exactly one of: 'roc_length', 'ovl', or 'joint_method'")
    }

    ## -----------------------------
    ## 7. Check "group" values
    ##    & filter data based on that
    ## -----------------------------
    if (!all(group %in% c(1, 2)) || !all(c(1, 2) %in% group)) {
      stop("'group' must contain exactly two groups coded as 1 and 2.")
    }

    ## -----------------------------
    ## 8. Check progress bar
    ## -----------------------------
    if (!is.logical(progress) || length(progress) != 1) {
      stop("'progress' must be either TRUE or FALSE.")
    }

    ## -----------------------------
    ## 9. Check plot argument
    ## -----------------------------
    if (!is.logical(plot) || length(plot) != 1) {
      stop("'plot' must be either TRUE or FALSE.")
    }

    if (plot && method != "joint_method") {
      stop("Plotting is only available for method = 'joint_method'. Set plot = FALSE or change method.")
    }

    time1 <- time[group == 1]
    censor1 <- status[group == 1]
    time2 <- time[group == 2]
    censor2 <- status[group == 2]

    ## -----------------------------
    ## Run selected method
    ## -----------------------------
    if (method == "roc_length") {
      return(roc_LeftCenSurvival_test(time1, censor1, time2, censor2, n_perm, progress))
    }

    if (method == "ovl") {
      return(ovl_LeftCenSurvival_test(time1, censor1, time2, censor2, n_perm, progress))
    }

    if (method == "joint_method") {
      return(joint.roc_ovl_LeftCenSurvival_test(time1, censor1, time2, censor2, n_perm, progress, plot))
    }
  }


  else if (censor_type == "double") {

    ## -----------------------------
    ## 1. Check lengths
    ## -----------------------------
    if (length(unique(c(length(time), length(status), length(group)))) != 1) {
      stop("Length mismatch: 'time', 'status', and 'group' must have the same length.")
    }

    ## -----------------------------
    ## 2. Check missing values
    ## -----------------------------
    if (any(is.na(c(time, status, group)))) {
      stop("Missing values found in at least one of 'time', 'status', or 'group'.")
    }

    ## -----------------------------
    ## 3. Check data types
    ## -----------------------------
    if (!is.numeric(time)) {
      stop("'time' must be a numeric vector.")
    }

    if (!is.numeric(status)) {
      stop("'status' must be a numeric vector.")
    }

    if (!is.numeric(group)) {
      stop("'group' must be a numeric vector.")
    }

    ## -----------------------------
    ## 4. Check censoring values
    ## -----------------------------
    if (!all(status %in% c(0, 1, -1)) || !all(c(0, 1, -1) %in% status)) {
      stop("'status' must contain all three values: 0 = event, 1 = right-censored, and -1 = left-censored, when censor_type = 'double'.")
    }

    ## -----------------------------
    ## 5. Check n_perm
    ## -----------------------------
    if (!is.numeric(n_perm) || length(n_perm) != 1 || n_perm <= 0 || n_perm %% 1 != 0) {
      stop("'n_perm' must be a positive integer.")
    }

    ## -----------------------------
    ## 6. Check method
    ## -----------------------------
    valid_methods <- c("roc_length", "ovl", "joint_method")

    if (length(method) != 1 || !method %in% valid_methods) {
      stop("'method' must be exactly one of: 'roc_length', 'ovl', or 'joint_method'")
    }

    ## -----------------------------
    ## 7. Check "group" values
    ##    & filter data based on that
    ## -----------------------------
    if (!all(group %in% c(1, 2)) || !all(c(1, 2) %in% group)) {
      stop("'group' must contain exactly two groups coded as 1 and 2.")
    }

    ## -----------------------------
    ## 8. Check progress bar
    ## -----------------------------
    if (!is.logical(progress) || length(progress) != 1) {
      stop("'progress' must be either TRUE or FALSE.")
    }

    ## -----------------------------
    ## 9. Check plot argument
    ## -----------------------------
    if (!is.logical(plot) || length(plot) != 1) {
      stop("'plot' must be either TRUE or FALSE.")
    }

    if (plot && method != "joint_method") {
      stop("Plotting is only available for method = 'joint_method'. Set plot = FALSE or change method.")
    }

    time1 <- time[group == 1]
    censor1 <- status[group == 1]
    time2 <- time[group == 2]
    censor2 <- status[group == 2]

    ## -----------------------------
    ## Run selected method
    ## -----------------------------
    if (method == "roc_length") {
      return(roc_DoubleCenSurvival_test(time1, censor1, time2, censor2, n_perm, progress))
    }

    if (method == "ovl") {
      return(ovl_DoubleCenSurvival_test(time1, censor1, time2, censor2, n_perm, progress))
    }

    if (method == "joint_method") {
      return(joint.roc_ovl_DoubleCenSurvival_test(time1, censor1, time2, censor2, n_perm, progress, plot))
    }
  }

  else
    stop("'censor_type' must be exactly one of: 'right', 'left', or 'double'")

}

