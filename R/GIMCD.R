#' Gaussian imputation followed by MCD
#'
#' Gaussian imputation uses the classical non-robust mean and covariance
#' estimator and then imputes predictions under the multivariate normal model.
#' Outliers may be created by this procedure. Then a high-breakdown robust
#' estimate of the location and scatter with the Minimum Covariance Determinant
#' algorithm is obtained and finally outliers are determined based on Mahalanobis
#' distances based on the robust location and scatter.
#'
#' Normal imputation from package \code{norm} and MCD from package \code{MASS}.
#' Note that currently MCD does not accept weights.
#'
#' @param data a data frame or matrix with the data.
#' @param alpha a threshold value for the cut-off for the outlier
#' Mahalanobis distances.
#' @param seedem random number generator seed for EM algorithm
#' @param seedmcd random number generator seed for MCD algorithm,
#' if \code{seedmcd} is missing, an internal seed will be used.
#' @return Result is stored in a global list GIMCD.r:
#' \describe{
#'   \item{\code{center}}{robust center}
#'   \item{\code{scatter}}{robust covariance}
#'   \item{\code{alpha}}{quantile for cut-off value}
#'   \item{\code{computation.time}}{elapsed computation time}
#'   \item{\code{outind}}{logical vector of outlier indicators}
#'   \item{\code{dist}}{Mahalanobis distances}
#' }
#' @author Beat Hulliger
#' @references Béguin, C. and Hulliger, B. (2008), The BACON-EEM Algorithm
#' for Multivariate Outlier Detection, in Incomplete Survey Data, Survey
#' Methodology, Vol. 34, No. 1, pp. 91-103.
#' @seealso \code{\link[MASS]{cov.rob}}
#' @examples
#' data(bushfirem)
#' det.res <- GIMCD(bushfirem, alpha = 0.1)
#' print(det.res$center)
#' PlotMD(det.res$dist, ncol(bushfirem))
#' @export
#' @importFrom stats mahalanobis median qf
GIMCD <- function(data, alpha = 0.05, seedem = 23456789, seedmcd) {

  # start computation time
  calc.time <- proc.time()

  # set seed
  norm::rngseed(seedem)

  # transform to matrix
  if (!is.matrix(data)) {data <- as.matrix(data)}

  # correct alpha if alpha smaller 0.5
  if (alpha < 0.5) {alpha <- 1 - alpha}

  # sorts rows of data by missingness patterns, and centers/scales columns
  s <- norm::prelim.norm(data)

  # get normal parameters with the help of EM algorithm
  thetahat <- norm::em.norm(s, showits = FALSE)

  # impute missing multivariate normal data
  imp.data <- norm::imp.norm(s, thetahat, data)

  ### MCD algorithm

  # set seed if not missing
  if (!missing(seedmcd)) {set.seed(seedmcd)}

  # computes multivariate location and scale estimate
  MCD.cov <- MASS::cov.mcd(imp.data)

  # returns the squared Mahalanobis distance
  dist <- mahalanobis(imp.data, MCD.cov$center, MCD.cov$cov)

  # number of rows
  n <- nrow(data)

  # number of columns
  p <- ncol(data)

  # compute cutpoint based on Mahalanobis distances
  cutpoint <- qf(alpha, p, n - p) * median(dist) / qf(0.5, p, n - p)

  # get indices of outliers
  outind <- (dist > cutpoint)

  # stop computation time
  calc.time <- proc.time() - calc.time

  # output to console
  cat("GIMCD has detected", sum(outind), "outliers in", calc.time[1], "seconds.")

  # return output
  return(list(
    center = MCD.cov$center,
    scatter = MCD.cov$cov,
    alpha = 1 - alpha,
    computation.time = calc.time[1],
    cutpoint = cutpoint,
    outind = outind,
    dist = dist))
}
