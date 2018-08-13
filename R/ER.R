#' Robust EM-algorithm ER
#'
#' The \code{ER} function is an implementation of the ER-algorithm
#' of Little and Smith (1987).
#'
#' The M-step of the EM-algorithm uses a one-step M-estimator.
#'
#' @param data a data frame or matrix with the data.
#' @param weights sampling weights.
#' @param alpha probability for the quantile of the cut-off.
#' @param psi.par further parameters passed to the psi-function.
#' @param em.steps number of iteration steps of the EM-algorithm.
#' @param steps.output if \code{TRUE}, verbose output.
#' @param Estep.output if \code{TRUE}, estimators are output at each iteration.
#' @param tolerance convergence criterion (relative change).
#' @return
#' \describe{
#'   \item{\code{sample.size}}{Number of observations}
#'   \item{\code{number.of.variables}}{Number of variables}
#'   \item{\code{significance.level}}{alpha}
#'   \item{\code{computation.time}}{Elapsed computation time}
#'   \item{\code{good.data}}{Indices of the data in the final good subset}
#'   \item{\code{outliers}}{Indices of the outliers}
#'   \item{\code{center}}{Final estimate of the center}
#'   \item{\code{scatter}}{Final estimate of the covariance matrix}
#'   \item{\code{dist}}{Final Mahalanobis distances}
#'   \item{\code{rob.weights}}{Robustness weights in the final EM step}
#' }
#' @author Beat Hulliger
#' @references Little, R. and P. Smith (1987). Editing and imputation for
#' quantitative survey data. Journal of the American Statistical Association, 82, 58-68.
#' @seealso \code{\link{BEM}}
#' @examples
#' data(bushfirem, bushfire.weights)
#' det.res <- ER(bushfirem, weights = bushfire.weights, alpha = 0.05,
#' steps.output = TRUE, em.steps = 100, tol = 2e-6)
#' PlotMD(det.res$dist, ncol(bushfirem))
#' @export
#' @importFrom stats cov.wt weighted.mean qchisq
ER <- function(data, weights, alpha = 0.01, psi.par = c(2, 1.25), em.steps = 100,
               steps.output = FALSE, Estep.output = FALSE, tolerance = 1e-6) {

  # ------- preprocessing of the data -------

  # Removing the unit(s) with all items missing

  # transform to matrix
  if (!is.matrix(data)) {data <- as.matrix(data)}

  # number of rows
  n <- nrow(data)

  # number of columns
  p <- ncol(data)

  # correct alpha if alpha smaller 0.5
  if (alpha < 0.5) {alpha <- 1 - alpha}

  # set weights to 1 if missing
  if (missing(weights)) {weights <- rep(1, n)}

  # get indices of observations where at least one column is not missing
  new.indices <- which(apply(is.na(data), 1, prod) == 0)

  # remove completely missing observations
  if (length(new.indices) < n) {
    cat("Warning: missing observations", which(apply(is.na(data), 1, prod) == 1),
        "removed from the data\n")
    data <- data[new.indices, ]
    weights <- weights[new.indices]
    n <- nrow(data)
  }

  # print progress to console
  if (steps.output) {cat("End of preprocessing\n")}

  # ------- start computation -------

  # start computation time
  calc.time <- proc.time()

  # create missingness patterns per row (e.g. "00100")
  s.patterns <- apply(matrix(as.integer(is.na(data)), n, p),
                      1, paste, sep = "", collapse = "")

  # order rows depending on pattern (starting with e.g. "00000")
  perm <- order(s.patterns)

  # resort data, patterns, and weights
  data <- data[perm, ]
  s.patterns <- s.patterns[perm]
  weights <- weights[perm]

  # count occurences per pattern
  s.counts <- as.vector(table(s.patterns))

  # indices of the last observation of each missingness pattern in the
  # dataset ordered by missingness pattern.
  s.id <- cumsum(s.counts)

  # total number of different missingness patterns
  S <- length(s.id)

  # missing items for each pattern
  missing.items <- is.na(data[s.id, , drop = FALSE])

  # number of missing items for each pattern
  nb.missing.items <- apply(missing.items, 1, sum)

  # print progress to console
  if (steps.output) {cat("End of missingness statistics\n")}

  # ------- preparation for call to ER -------

  # Case where some observations are complete (no missing item) =>
  # start.mean and start.var are calculated on the complete observations
  if (nb.missing.items[1] == 0) {

    # returns weighted covariance matrix and mean
    cov.complete <- cov.wt(data[1:s.counts[1], ],
                           wt = weights[1:s.counts[1]] / sum(weights[1:s.counts[1]]))

    # get mean and cov
    mean.start <- cov.complete$center
    var.start <- cov.complete$cov

  } else {
    # Case where all observations have missing items
    mean.start <- apply(data, 2, weighted.mean, w = weights, na.rm = TRUE)
    var.start <- diag(apply(data, 2, weighted.var, w = weights, na.rm = TRUE))
  }

  # ------- mean and covariance matrix computed by ER -------

  # print progress to console
  if (steps.output) {
    cat("\n", "start.mean: ", mean.start, "\n", "start.var:\n")
    print(var.start)
    cat("\n")
  }

  # call internal function to compute mean and covariance matrix
  ER.result <- .ER.normal(data = data, weights, psi.par = psi.par, np = sum(weights),
                          p = p, s.counts = s.counts, s.id = s.id, S,
                          missing.items = missing.items, nb.missing.items,
                          start.mean = mean.start, start.var = var.start,
                          numb.it = em.steps, Estep.output = Estep.output,
                          tolerance = tolerance)

  # get mean and cov
  ER.mean <- ER.result$theta[1, 2:(p + 1)]
  ER.var <- ER.result$theta[2:(p + 1), 2:(p + 1)]

  # stop computation time
  calc.time <- proc.time() - calc.time

  # Nominate the outliers using the original numbering
  dist <- ER.result$dist
  good <- dist <= qchisq(alpha, p)
  outliers <- perm[!good]

  # output to console
  cat("\n", "ER has detected", sum(!good), "outlier(s) in", calc.time[3], "seconds.", "\n", "\n")
  if (!ER.result$convergence) {cat("\n", "ER did not converge.", "\n")}

  # return output
  return(list(
    sample.size = n,
    number.of.variables = p,
    significance.level = alpha,
    computation.time = calc.time,
    good.data = perm[good],
    outliers = outliers,
    center = ER.mean,
    scatter = ER.var,
    dist = dist[order(perm)],
    rob.weights = ER.result$rob.weights))
}

