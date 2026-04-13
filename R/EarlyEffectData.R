#' Early Effect Survival Dataset (Gastrointestinal Tumor Study Group)
#'
#' A survival dataset derived from a randomized clinical trial conducted by the
#' Gastrointestinal Tumor Study Group, comparing two treatment strategies for
#' patients with locally advanced gastric cancer.
#'
#' The primary objective of the study was to compare survival outcomes between:
#' \itemize{
#'   \item Chemotherapy alone
#'   \item Chemotherapy combined with radiation therapy
#' }
#'
#' Kaplan–Meier survival curves initially show higher survival probabilities for
#' patients receiving chemotherapy alone. However, the hazard in this group increases
#' over time, and the survival curves cross at later time points, indicating a
#' non-proportional hazards scenario. This pattern represents a classic example of
#' an early treatment effect with diminishing long-term differences.
#'
#' @format A data frame with 90 observations and 3 variables:
#' \describe{
#'   \item{time}{Observed survival time}
#'   \item{status}{Event indicator (0 = event occurred, 1 = right-censored)}
#'   \item{group}{Treatment group (1 = chemotherapy alone, 2 = combined therapy)}
#' }
#'
#' @details
#' The dataset contains two groups of equal size (n = 45 each). It is widely used
#' as an example of survival data with crossing hazards, where the proportional
#' hazards assumption is violated.
#'
#' Note that the censoring indicator follows the convention:
#' \itemize{
#'   \item 0 = event (failure)
#'   \item 1 = right censoring
#' }
#'
#' @source
#' Stablein DM, Carter WH Jr, Novak JW (1981).
#' Analysis of survival data with nonproportional hazard functions.
#' \emph{Controlled Clinical Trials}, 2(2), 149–159.
"EarlyEffectData"
