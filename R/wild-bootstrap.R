#----------------------------------------------
# auxilliary distribution functions
#----------------------------------------------

r_Mammen <- function(n) {
  pts <- c(-(sqrt(5) - 1), sqrt(5) + 1) / 2
  prob <- (sqrt(5) + 1) / (2 * sqrt(5))
  sample(x = pts, size = n, replace = TRUE, prob = c(prob, 1 - prob))
}

r_Rademacher <- function(n) {
  sample(x = c(-1L, 1L), size = n, replace = TRUE)
}

r_sixpoint <- function(n) {
  six_points <- c(-sqrt(3/2), -1, -sqrt(1 / 2), sqrt(1 / 2), 1, sqrt(3 / 2))
  sample(x = six_points, size = n, replace = TRUE)
}

#----------------------------------------------
# wild bootstrap function
#----------------------------------------------

#' Cluster-wild bootstrap F test
#'
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix
#' of a set of regression coefficient estimates.
#'
#' @param obj Fitted null model from which to calculate residuals for use in the
#'   cluster-wild bootstrap algorithm.
#' @param constraints formula specifying the set of additional parameter
#'   restrictions to be tested.
#' @param bootstraps Integer specifying number of bootstrap replicates to
#'   generate.
#' @param auxilliary_dist Function for generating auxilliary random variables
#'   for perturbing the residuals in the cluster wild bootstrap algorithm. The
#'   default is \code{r_Rademacher}, which samples from a two-point distribution
#'   with equal mass at -1 and 1. A user-specified function may be passed
#'   instead.
#' @param target Optional matrix or vector describing the working
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4}
#'   adjustment matrices. If a vector, the target matrix is assumed to be
#'   diagonal. If not specified, \code{vcovCR} will attempt to infer a value.
#' @residual_adjustment Character string specifying which small-sample
#'   adjustment to use when estimating the residuals under the null model.
#'   Available options are the same as the \code{type} argument of
#'   \code{\link{vcovCR}}. Default is \code{"CR2"}.
#' @test_adjustment. Character string specifying which small-sample adjustment
#' to use when calculating the Wald test statistic. Available options are the
#' same as the \code{type} argument of \code{\link{vcovCR}}. Default is
#' \code{"CR0"}.
#' @param ... Additional arguments passed to \code{\link{vcovCR}} when
#'   calculating the adjustment matrices used to estimate the residauls.
#' @inheritParams vcovCR
#'
#' @description
#' @details
#' @references
#' @return
#' @examples
#'
#'
#' @export


Wild_bootstrap <- function(obj, constraints,
                           cluster, bootstraps = 2000,
                           auxilliary_dist = r_Rademacher,
                           residual_adjustment = "CR2",
                           test_adjustment = "CR0",
                           ...) {

  # fit expanded model by updating obj according to constraints

  # calculate Wald test for updated model
  vcov_updated <-
  Qstat <- as.numeric(t(beta) %*% chol2inv(chol(vcov_updated)) %*% beta)

  # prepare for bootstrapping

  J <- nlevels(cluster)
  vcov_null <- vcovCR(obj, cluster = cluster, type = residual_adjustment, ...)
  res_list_null <- split(residuals_CS(obj), cluster)
  B_mats <- attr(vcov_null, "adjustments")
  eta_j <- auxilliary_dist(n = J)

  # simulate bootstrap distribution of Wald test statistic

  # calculate p-value from bootstrapped distribution

  res <- data.frame(Wald_statistic = Q, p_val = p_val)
  return(res)
}
