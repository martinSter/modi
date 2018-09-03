#' BACON-EEM Algorithm for multivariate outlier detection in incomplete
#' multivariate survey data
#'
#' \code{BEM} starts from a set of uncontaminated data with possible
#' missing values, applies a version of the EM-algorithm to estimate
#' the center and scatter of the good data, then adds (or deletes)
#' observations to the good data which have a Mahalanobis distance
#' below a threshold. This process iterates until the good data remain
#' stable. Observations not among the good data are outliers.
#'
#' The BACON algorithm with \code{v = 1} is not robust but affine equivariant
#' while \code{v = 1} is robust but not affine equivariant. The threshold for
#' the (squared) Mahalanobis distances, beyond which an observation is an
#' outlier, is a standardised chisquare quantile at \code{(1 - alpha)}. For
#' large data sets it may be better to choose \code{alpha / n} instead. The
#' internal function \code{.EM.normal} is usually called from \code{BEM}.
#' \code{.EM.normal} is implementing the EM-algorithm in such a way that
#' part of the calculations can be saved to be reused in the \code{BEM}
#' algorithm. \code{.EM.normal} does not contain the computation of the
#' observed sufficient statistics, they will be computed in the main
#' program of \code{BEM} and passed as parameters as well as the statistics
#' on the missingness patterns.
#'
#' @param data a matrix or data frame. As usual, rows are observations and
#' columns are variables.
#' @param weights a non-negative and non-zero vector of weights for each
#' observation. Its length must equal the number of rows of the data.
#' Default is \code{rep(1, nrow(data))}.
#' @param v an integer indicating the distance for the definition of the
#' starting good subset: \code{v = 1} uses the Mahalanobis distance based
#' on the weighted mean and covariance, \code{v = 2} uses the Euclidean
#' distance from the componentwise median.
#' @param c0 the size of initial subset is \code{c0 * ncol(data)}.
#' @param alpha a small probability indicating the level \code{(1 - alpha)}
#' of the cutoff quantile for good observations.
#' @param md.type type of Mahalanobis distance: \code{"m"} marginal,
#' \code{"c"} conditional.
#' @param em.steps.start number of iterations of EM-algorithm for starting
#' good subset.
#' @param em.steps.loop number of iterations of EM-algorithm for good subset.
#' @param better.estimation if \code{better.estimation = TRUE}, then the
#' EM-algorithm for the final good subset iterates \code{em.steps.start} more.
#' @param monitor if \code{TRUE}, verbose output.
#' @return \code{BEM} returns a list whose first component \code{output} is a
#' sublist with the following components:
#' \describe{
#'   \item{\code{sample.size}}{Number of observations}
#'   \item{\code{discarded.observations}}{Number of discarded observations}
#'   \item{\code{number.of.variables}}{Number of variables}
#'   \item{\code{significance.level}}{The probability used for the cutpoint,
#'   i.e. \code{alpha}}
#'   \item{\code{initial.basic.subset.size}}{Size of initial good subset}
#'   \item{\code{final.basic.subset.size}}{Size of final good subset}
#'   \item{\code{number.of.iterations}}{Number of iterations of the BACON step}
#'   \item{\code{computation.time}}{Elapsed computation time}
#'   \item{\code{center}}{Final estimate of the center}
#'   \item{\code{scatter}}{Final estimate of the covariance matrix}
#'   \item{\code{cutpoint}}{The threshold MD-value for the cut-off of outliers}
#' }
#' The further components returned by \code{BEM} are:
#' \describe{
#'   \item{\code{outind}}{Indicator of outliers}
#'   \item{\code{dist}}{Final Mahalanobis distances}
#' }
#' @author Beat Hulliger
#' @note \code{BEM} uses an adapted version of the EM-algorithm in function
#' \code{.EM-normal}.
#' @references Béguin, C. and Hulliger, B. (2008) The BACON-EEM Algorithm for
#' Multivariate Outlier Detection in Incomplete Survey Data, Survey Methodology,
#' Vol. 34, No. 1, pp. 91-103.
#'
#' Billor, N., Hadi, A.S. and Vellemann, P.F. (2000). BACON: Blocked Adaptative
#' Computationally-efficient Outlier Nominators. Computational Statistics and
#' Data Analysis, 34(3), 279-298.
#'
#' Schafer J.L. (2000), Analysis of Incomplete Multivariate Data, Monographs on
#' Statistics and Applied Probability 72, Chapman & Hall.
#' @examples
#' # Bushfire data set with 20% MCAR
#' data(bushfirem, bushfire.weights)
#' bem.res <- BEM(bushfirem, bushfire.weights, alpha = (1 - 0.01 / nrow(bushfirem)))
#' print(bem.res$output)
#' @export
#' @importFrom stats qchisq cov.wt mahalanobis
BEM <- function(data, weights, v = 2, c0 = 3, alpha = 0.01, md.type = "m",
                em.steps.start = 10, em.steps.loop = 5,
                better.estimation = FALSE, monitor = FALSE) {

  # ------- preprocessing -------

  # transform data to matrix
  if(!is.matrix(data)) {data <- as.matrix(data)}

  # number of rows
  n <- nrow(data)

  # number of columns
  p <- ncol(data)

  # set weights to 1 if missing
  if (missing(weights)) {(weights <- rep(1, n))}

  # finding the unit(s) with all items missing
  new.indices <- which(apply(is.na(data), 1, prod) == 0)
  discarded <- NA
  nfull <- n

  # remove observations with all missings
  if (length(new.indices) < n) {
    discarded <- which(apply(is.na(data), 1, prod) == 1)
    cat("Warning: missing observations", discarded, "removed from the data\n")
    data <- data[new.indices, ]
    weights <- weights[new.indices]
    n <- nrow(data)
  }

  # print progress to console
  if (monitor) {cat("End of preprocessing\n")}

  # start computation time
  calc.time <- proc.time()

  # Order the data by missingness patterns :
  # s.patterns = vector of length n, with the
  # missingness patterns stocked as strings of
  # the type "11010...011" with "1" for missing.
  # data, weights and s.patterns ordered using
  # s.patterns' order
  s.patterns <- apply(matrix(as.integer(is.na(data)), n, p),
                      1, paste, sep = "", collapse = "")

  perm <- order(s.patterns)
  data <- data[perm, ]
  s.patterns <- s.patterns[perm]
  weights <- weights[perm]

  # Missingness patterns stats :
  # s.counts = counts of the different missingness
  # patterns ordered alphabetically.
  # s.id = indices of the last observation of each
  # missingness pattern in the dataset ordered by
  # missingness pattern.
  # S = total number of different missingness patterns
  # missing.items = missing items for each pattern
  # nb.missing.items = number of missing items for each pattern
  s.counts <- as.vector(table(s.patterns))
  s.id <- cumsum(s.counts)
  S <- length(s.id)
  missing.items <- is.na(data[s.id, , drop = FALSE])
  nb.missing.items <- apply(missing.items, 1, sum)

  # print progress to console
  if (monitor) {cat("End of missingness statistics\n")}

  # ------- constants -------
  # constants used by the BACON algorithm
  # to select the good points

  # ???
  N <- sum(weights)
  initial.length <- c0 * p

  # stop program if inital.length too large
  if (initial.length > n) {
    stop("\nInitial length bigger than number of
         observations. Please, decrease c0.")
  }

  # ???
  c.np <- 1 + (p + 1) / (N - p) + 2 / (N - 1 - 3 * p)
  h <- floor((N + p + 1) / 2)

  # correct alpha if alpha smaller 0.5
  if (alpha > 0.5) {
    alpha <- 1 - alpha
    cat("alpha should be less than 0.5:
        alpha set to 1 - alpha\n", 1 - alpha)
  }

  # determine 'chi.sq'
  chi.sq <- qchisq(1 - alpha, p)

  # ------- STEP 1 -------
  # The two possible starts for BACON, modified
  # to deal with missing items: Version 2 (default):
  # the distances are the Euclidean distances form the
  # componentwise median; the median is computed by
  # removing the missing items in each variable;
  # missing values in an observation are not included
  # in the distance to the median. If pg is the number
  # of columns in which no missing values occur for that
  # observation, then the distance returned is sqrt(p/pg)
  # times the Euclidean distance between the two vectors
  # of length pg shortened to exclude missing values.
  # Version 1 : the usual mean and covariance matrix are
  # used to compute Mahalanobis distances; both mean and
  # covariance are computed by EM if any missing value
  # occurs; the Mahalanobis distance for an observation is
  # computed by restricting the mean and the covariance
  # matrix to the subspace of non-missing variables and
  # the sqrt(p/pg) correction is used as in Version 1.

  if (v == 2) {

    # version 2
    EM.mean <- apply(data, 2, weighted.quantile, w = weights)
    dist <- apply(sweep(data, 2, EM.mean) ^ 2, 1, sum, na.rm = TRUE) *
      p / (p - apply(is.na(data), 1, sum))

    } else {

      # version 1
      if (S == 1 & nb.missing.items[1] == 0) {

        # case where no missing value occurs =>
        # regular mean and covariance matrix
        EM.mean <- apply(data, 2, weighted.mean, w = weights)
        EM.var <- cov.wt(data, wt = weights, center = EM.mean, method = "ML")$cov

        } else {

          # case where missing values occur =>
          # mean and covariance matrix computed by EM
          T.obs <- matrix(0, p + 1, p + 1)

          if (nb.missing.items[1] == 0) {

            # case where some observations are complete
            # (no missing item) => computation of the
            # sufficient statistics on these observations
            # and then EM.

            # print progress to console
            if (monitor) {cat("Preparation of T_obs for version 1\n")}

            weights.obs <- weights[1:s.counts[1]]
            T.obs[1, ] <- T.obs[ ,1] <-
              c(sum(weights.obs), apply(weights.obs * data[1:s.counts[1], ], 2, sum))

            for (i in 1:s.counts[1]) {

              T.obs[2:(p + 1), 2:(p + 1)] <-
                T.obs[2:(p + 1), 2:(p + 1)] + weights.obs[i] * data[i, ] %*% t(data[i, ])

            }

            EM.result <- .EM.normal(data = data[(s.id[1] + 1):n, , drop = FALSE],
                                    weights = weights[(s.id[1] + 1):n], n = N, p = p,
                                    s.counts = s.counts[2:S], s.id = s.id[2:S] - s.id[1],
                                    S = S - 1, T.obs = T.obs,
                                    start.mean = apply(data, 2, weighted.mean,
                                                       w = weights, na.rm = TRUE),
                                    start.var = diag(apply(data, 2, weighted.var,
                                                           w = weights, na.rm = TRUE)),
                                    numb.it = em.steps.start, Estep.output = monitor)

            } else {

              # case where all observations have missing items => EM
              EM.result <- .EM.normal(data = data, weights, n = N, p = p,
                                      s.counts = s.counts, s.id = s.id, S,
                                      T.obs = T.obs,
                                      start.mean = apply(data, 2, weighted.mean,
                                                         w = weights, na.rm = TRUE),
                                      start.var = diag(apply(data, 2, weighted.var,
                                                             w = weights, na.rm = TRUE)),
                                      numb.it = em.steps.start, Estep.output = monitor)

            }

          EM.mean <- EM.result[1, 2:(p + 1)]
          EM.var <- EM.result[2:(p + 1), 2:(p + 1)]

        }

      # prepare indices (used in comp. of Mahalanobis dist.)
      indices <- (!missing.items[1, ])

      # type of Mahalanobis distances (c = conditional)
      if (md.type == "c") {
        EM.var.inverse <- solve(EM.var)
        } else {
          EM.var.inverse <- EM.var
        }

      # computation of the Mahalanobis distances
      dist <- mahalanobis(data[1:s.id[1], indices, drop = FALSE],
                          EM.mean[indices], EM.var.inverse[indices, indices],
                          inverted = (md.type == "c")) * p / (p - nb.missing.items[1])

      # Mahalanobis distances if S > 1
      if (S > 1) {

        for (i in 2:S) {

          indices <- (!missing.items[i, ])
          dist <- c(dist,
                    mahalanobis(data[(s.id[i-1]+1):s.id[i], indices, drop = FALSE],
                                EM.mean[indices],
                                EM.var.inverse[indices, indices, drop = FALSE],
                                inverted = (md.type == "c")) * p / (p - nb.missing.items[i]))

        }

      }

    }

  # print progress to console
  if (monitor) {cat("Version ", v, ": estimation of mean = ", signif(EM.mean, 4), "\n")}

  # selection of the initial basic good subset using the computed distance
  ordre <- order(dist)
  good <- ordre[1:initial.length]
  good <- good[order(good)]

  # print progress to console
  if (monitor) {cat("Initial good subset size: ", length(good), "\n")}

  # statistics of the missingness patterns of the good subset
  n.good <- length(good)
  s.patterns.good <- s.patterns[good]
  s.counts.good <- as.vector(table(s.patterns.good))
  s.id.good <- cumsum(s.counts.good)
  S.good <- length(s.id.good)
  missing.items.good <- is.na(data[good[s.id.good], , drop = FALSE])
  nb.missing.items.good <- apply(missing.items.good, 1, sum)
  weights.good <- weights[good]
  N.good <- sum(weights.good)

  # determination of the indices (if any) of the observations
  # in the good subset without missing items
  T.obs.good.exist <- (nb.missing.items.good[1] == 0)

  if (T.obs.good.exist) {
    T.obs.good.indices <- good[1:s.id.good[1]]
  }

  # computation of mean and covariance matrix of the good subset by EM
  if (S.good == 1 & T.obs.good.exist) {

    # case where no missing value occurs => regular mean and
    # covariance matrix computed and the sufficient statistics
    # deduced from them
    EM.mean.good <- apply(data[good, ], 2, weighted.mean, w = weights.good)
    EM.var.good <-
      cov.wt(data[good, ], wt = weights.good, center = EM.mean.good, method = "ML")$cov

    T.obs.good <- matrix(0, p + 1, p + 1)
    T.obs.good[1, ] <- T.obs.good[, 1] <- c(-1, EM.mean.good)
    T.obs.good[2:(p + 1), 2:(p + 1)] <- EM.var.good * (N.good - 1) / N.good
    T.obs.good <- .sweep.operator(T.obs.good, 1, TRUE) * N.good
    T.obs.good.exist <- TRUE

    } else {

      # case where missing values occur =>
      # mean and covariance matrix computed by EM
      T.obs.good <- matrix(0, p + 1, p + 1)

      if (T.obs.good.exist) {

        # case where some observations are complete (no missing item) =>
        # computation of the sufficient statistics on these observations
        # and then EM.

        # print progress to console
        if (monitor) {cat("Preparation of T_obs\n")}

        weights.good.obs <- weights.good[1:s.counts.good[1]]
        T.obs.good[1, ] <- T.obs.good[, 1] <-
          c(sum(weights.good.obs),
            apply(weights.good.obs * data[good[1:s.counts.good[1]], , drop = FALSE], 2, sum))

        # start loop
        for (i in 1:s.counts.good[1]) {

          T.obs.good[2:(p + 1), 2:(p + 1)] <- T.obs.good[2:(p + 1), 2:(p + 1)] +
            weights.good.obs[i] * data[good[i], ] %*% t(data[good[i], ])

        }

        EM.result.good <-
          .EM.normal(data = data[good[(s.id.good[1] + 1):n.good], , drop = FALSE],
                     weights = weights.good[(s.id.good[1] + 1):n.good], n = N.good, p = p,
                     s.counts = s.counts.good[2:S.good],
                     s.id = s.id.good[2:S.good] - s.id.good[1],
                     S = S.good - 1, T.obs = T.obs.good,
                     start.mean = apply(data[good, ], 2, weighted.mean,
                                        w = weights.good, na.rm = TRUE),
                     start.var = diag(apply(data[good, ], 2, weighted.var,
                                            w = weights.good, na.rm = TRUE)),
                     numb.it = em.steps.loop, Estep.output = monitor)

        } else {

          # case where all observations have missing items => EM
          EM.result.good <-
            .EM.normal(data = data[good, , drop = FALSE],
                       weights = weights.good, n = N.good, p = p,
                       s.counts = s.counts.good, s.id = s.id.good,
                       S = S.good, T.obs = T.obs.good,
                       start.mean = apply(data[good, ], 2, weighted.mean,
                                          w = weights.good, na.rm = TRUE),
                       start.var = diag(apply(data[good, ], 2, weighted.var,
                                              w = weights.good, na.rm = TRUE)),
                       numb.it = em.steps.loop, Estep.output = monitor)

        }

      # center and variance
      EM.mean.good <- EM.result.good[1, 2:(p + 1)]
      EM.var.good <- EM.result.good[2:(p + 1), 2:(p + 1)]

    }

  # print progress to console
  if (monitor) {
    cat("First good subset estimation of mean = ", signif(EM.mean.good, 4), "\n")
  }

  # test if the size of the good subset is not
  # too small (singular covariance matrix), otherwise stop program
  if (qr(EM.var.good)$rank < p) {
    stop("Initial subset size too small, please increase c0")
  }

  # ------- STEP 2 TO 4 -------
  # main loop of BACON: increase the good subset using
  # the Mahalanobis distances computed from it; the
  # Mahalanobis distances are computed as in Version 1
  # of the start (see above).

  # initialize 'count' variable
  count <- 0

  # start repeat loop
  repeat {

    # increase count
    count <- count + 1

    # upgrade of the constants values used by BACON
    r <- sum(weights.good)
    c.hr <- max(0, (h - r) / (h + r))
    test <- (c.np + c.hr) ^ 2 * chi.sq

    # prepare indices (used in comp. of Mahalanobis dist.)
    indices <- (!missing.items[1, ])

    # type of Mahalanobis distances (c = conditional)
    if (md.type == "c") {
      EM.var.good.inverse <- solve(EM.var.good)
    } else {
      EM.var.good.inverse <- EM.var.good
    }

    # computation of the Mahalanobis distances
    dist <- mahalanobis(data[1:s.id[1], indices, drop = FALSE],
                        EM.mean.good[indices],
                        EM.var.good.inverse[indices, indices],
                        inverted = (md.type == "c")) * p / (p - nb.missing.items[1])

    # Mahalanobis distances if S > 1
    if (S > 1) {

      for (i in 2:S) {

        indices <- (!missing.items[i, ])
        dist <- c(dist,
                  mahalanobis(data[(s.id[i - 1] + 1):s.id[i], indices, drop = FALSE],
                              EM.mean.good[indices],
                              EM.var.good.inverse[indices, indices, drop = FALSE],
                              inverted = (md.type == "c")) * p / (p - nb.missing.items[i]))

      }

    }

    # memorization of some data on the preceeding good subset
    oldgood <- good
    T.obs.oldgood.exist <- T.obs.good.exist

    if (T.obs.oldgood.exist) {
      T.obs.oldgood.indices <- T.obs.good.indices
    }

    # determination of the new good subset
    good <- (1:n)[dist <= test]

    # print progress to console
    if (monitor) {
      cat("Loop ", count, ": start; test = ", test, "; good subset =", length(good), "\n")
    }

    # comparaison with the preceeding good subset, break if equality
    if (length(good) == length(oldgood)) {
      if (prod(good == oldgood)) {
        break
      }
    }

    # statistics of the missingness patterns of the new good subset
    n.good <- length(good)
    s.patterns.good <- s.patterns[good]
    s.counts.good <- as.vector(table(s.patterns.good))
    s.id.good <- cumsum(s.counts.good)
    S.good <- length(s.id.good)
    missing.items.good <- is.na(data[good[s.id.good], , drop = FALSE])
    nb.missing.items.good <- apply(missing.items.good, 1, sum)
    weights.good <- weights[good]
    N.good <- sum(weights.good)

    # determination of the indices (if any) of the observations
    # in the new good subset without missing items
    T.obs.good.exist <- (nb.missing.items.good[1] == 0)

    if (T.obs.good.exist) {
      T.obs.good.indices <- good[1:s.id.good[1]]
    }

    # computation of mean and covariance matrix of the new good subset
    if (S.good == 1 & T.obs.good.exist) {

      # case where no missing value occurs =>
      # regular mean and covariance matrix computed
      # and the sufficient statistics deduced from them
      EM.mean.good <- apply(data[good, ], 2, weighted.mean, w = weights.good)
      EM.var.good <-
        cov.wt(data[good, ], wt = weights.good, center = EM.mean.good, method = "ML")$cov

      T.obs.good <- matrix(0, p + 1, p + 1)
      T.obs.good[1, ] <- T.obs.good[, 1] <- c(-1, EM.mean.good)
      T.obs.good[2:(p + 1), 2:(p + 1)] <- EM.var.good * (N.good - 1) / N.good
      T.obs.good <- .sweep.operator(T.obs.good, 1, TRUE) * N.good
      T.obs.good.exist <- TRUE

    } else {

      # case where missing values occur =>
      # mean and covariance matrix computed by EM with
      # starting parameters set to the preceeding estimates
      if (T.obs.good.exist) {

        # case where some observations are complete (no missing item) =>
        # computation of the sufficient statistics on these observations
        # and then EM with starting parameters set to the preceeding estimates
        if (T.obs.oldgood.exist) {

          # case where sufficient statistics from complete data were
          # computed in the preceeding good subset => upgrade of these
          # statistics substracting the observations that were deleted
          # from the precedding good subset and adding the observations
          # that were added to it

          # print progress to console
          if (monitor) {cat("Updating of T_obs\n")}

          good.boo <- oldgood.boo <- rep(FALSE, n)
          good.boo[T.obs.good.indices] <- TRUE
          oldgood.boo[T.obs.oldgood.indices] <- TRUE

          # start for loop
          for (i in (1:n)[xor(good.boo, oldgood.boo) & oldgood.boo]) {

            T.obs.good[1, ] <- T.obs.good[1, ] - weights[i] * c(1, data[i, ])
            T.obs.good[2:(p + 1), 2:(p + 1)] <-
              T.obs.good[2:(p + 1), 2:(p + 1)] - weights[i] * data[i, ] %*% t(data[i, ])

          }

          # start for loop
          for (i in (1:n)[xor(good.boo, oldgood.boo)&good.boo]) {

            T.obs.good[1, ] <- T.obs.good[1, ] + weights[i] * c(1, data[i, ])
            T.obs.good[2:(p + 1), 2:(p + 1)] <-
              T.obs.good[2:(p + 1), 2:(p + 1)] + weights[i] * data[i, ] %*% t(data[i, ])

          }

          T.obs.good[ ,1] <- T.obs.good[1, ]

          } else {

            # case where no sufficient statistics from complete
            # data were computed in the preceeding good subset =>
            # computation of these statistics for the new good subset

            # print progress to console
            if (monitor) {cat("New preparation of T_obs\n")}

            T.obs.good <- matrix(0, p + 1, p + 1)
            weights.good.obs <- weights.good[1:s.counts.good[1]]
            T.obs.good[1, ] <- T.obs.good[ ,1] <-
              c(sum(weights.good.obs),
                apply(weights.good.obs * data[good[1:s.counts.good[1]], , drop = F], 2,sum))

            # start for loop
            for (i in 1:s.counts.good[1]) {

              T.obs.good[2:(p + 1), 2:(p + 1)] <- T.obs.good[2:(p + 1), 2:(p + 1)] +
                weights.good.obs[i] * data[good[i], ] %*% t(data[good[i], ])

            }

          }

        EM.result.good <-
          .EM.normal(data = data[good[(s.id.good[1] + 1):n.good], , drop = F],
                     weights.good[(s.id.good[1] + 1):n.good], n = N.good, p = p,
                     s.counts = s.counts.good[2:S.good],
                     s.id = s.id.good[2:S.good] - s.id.good[1], S = S.good - 1,
                     T.obs = T.obs.good, start.mean = EM.mean.good,
                     start.var = EM.var.good, numb.it = em.steps.loop,
                     Estep.output = monitor)

      } else {

        # case where all observations have missing items =>
        # EM with starting parameters set to the preceeding estimates
        T.obs.good <- matrix(0, p + 1, p + 1)

        EM.result.good <-
          .EM.normal(data = data[good, , drop = F], weights.good,
                     n = N.good, p = p, s.counts = s.counts.good,
                     s.id = s.id.good, S = S.good, T.obs = T.obs.good,
                     start.mean = EM.mean.good, start.var = EM.var.good,
                     numb.it = em.steps.loop, Estep.output = monitor)

      }

      # center and variance
      EM.mean.good <- EM.result.good[1, 2:(p + 1)]
      EM.var.good <- EM.result.good[2:(p + 1), 2:(p + 1)]

    }

    # print progress to console
    if (monitor) {
      cat("Loop ", count, " end: estimation of mean = ", signif(EM.mean.good, 4), "\n")
    }

    # check if the computed convariance is singular, stop in that case.
    if (qr(EM.var.good)$rank < p) {
      stop("Singular covariance matrix with a particular subset (try a bigger alpha).")
    }

    # start next iteration of loop
    next

  }

  # stop computation time
  calc.time <- proc.time() - calc.time

  # ------- STEP 5 -------
  # nominate the outliers using the original
  # numbering (without discarded obs.)

  # determine outliers
  outliers <- perm[(1:n)[-good]]

  # if a better estimation is seeked, "em.steps.start"
  # more steps of EM are taken
  if (better.estimation) {

    if (T.obs.good.exist) {

      # case where some observations are complete,
      # more steps of EM with the sufficient
      # statistics computed before
      EM.result.good <- .EM.normal(data = data[good[(s.id.good[1] + 1):n.good], , drop = F],
                                   weights.good[(s.id.good[1] + 1):n.good], n = N.good,
                                   p = p, s.counts = s.counts.good[2:S.good],
                                   s.id = s.id.good[2:S.good] - s.id.good[1],
                                   S = S.good - 1, T.obs = T.obs.good,
                                   start.mean = EM.mean.good, start.var = EM.var.good,
                                   numb.it = em.steps.start, Estep.output = monitor)

      } else {

        # case where all observations have missing items =>
        # more steps of EM
        EM.result.good <- .EM.normal(data = data[good, , drop = F], weights.good,
                                     n = N.good, p = p, s.counts = s.counts.good,
                                     s.id = s.id.good, S = S.good, T.obs = T.obs.good,
                                     start.mean = EM.mean.good, start.var = EM.var.good,
                                     numb.it = em.steps.start, Estep.output = monitor)

      }

    EM.mean.good <- EM.result.good[1, 2:(p + 1)]
    EM.var.good <- EM.result.good[2:(p + 1), 2:(p + 1)]

    # prepare indices (used in comp. of Mahalanobis dist.)
    indices <- (!missing.items[1, ])

    # type of Mahalanobis distances (c = conditional)
    if (md.type == "c") {

      EM.var.good.inverse <- solve(EM.var.good)

      } else {

        EM.var.good.inverse <- EM.var.good

      }

    # computation of the Mahalanobis distances
    dist <- mahalanobis(data[1:s.id[1], indices, drop = F], EM.mean.good[indices],
                        EM.var.good.inverse[indices, indices],
                        inverted = (md.type == "c")) * p / (p - nb.missing.items[1])

    # Mahalanobis distances if S > 1
    if (S > 1) {

      for (i in 2:S) {

        indices <- (!missing.items[i, ])

        dist <-
          c(dist, mahalanobis(data[(s.id[i - 1] + 1):s.id[i], indices, drop = F],
                              EM.mean.good[indices],
                              EM.var.good.inverse[indices, indices, drop = F],
                              inverted = (md.type == "c")) * p / (p - nb.missing.items[i]))

      }

    }

  }

  # distances in the original numbering
  dist <- dist[order(perm)]
  cutpoint <- min(dist[outliers])

  # outliers and distances in original numbering with full dataset
  outn <- logical(n)
  outn[outliers] <- TRUE
  outnfull <- logical(nfull)
  outnfull[new.indices] <- outn
  distnfull <- rep(NA, nfull)
  distnfull[new.indices] <- dist

  # ------- results -------

  # prepare output
  BEM.r <- list(
    sample.size = n,
    discarded.observations = discarded,
    number.of.variables = p,
    significance.level = alpha,
    initial.basic.subset.size = initial.length,
    final.basic.subset.size = length(good),
    number.of.iterations = count,
    computation.time = calc.time,
    center = EM.mean.good,
    scatter = EM.var.good,
    cutpoint = cutpoint)

  BEM.i <- list(outind = outnfull, dist = distnfull)

	# output to console
	cat("\n", "BEM has detected", length(outliers),
	    "outlier(s) in", calc.time[1], "seconds.", "\n", "\n")

	# return output
	return(invisible(list(output = BEM.r, outind = outnfull,
	                      dist = distnfull)))
}

