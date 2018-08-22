#' Nearest Neighbour Imputation with Mahalanobis distance
#'
#' POEM takes into account missing values, outlier indicators, error indicators
#' and sampling weights.
#'
#' \code{POEM} assumes that an multivariate outlier detection has been carried out
#' beforehand and assumes the result is summarized in the vector \code{outind}.
#' In addition, further observations may have been flagged as failing edit-rules
#' and this information is given in the vector \code{errors}. The mean and
#' covariance estimate is calculated with the good observations (no outliers and
#' downweighted errors). Preliminary mean imputation is sometimes needed to avoid
#' a non-positive definite covariance estimate at this stage. Preliminary mean
#' imputation assumes that the problematic values of an observation (with errors,
#' outliers or missing) can be replaced by the mean of the rest of the non-problematic
#' observations. Note that the algorithm imputes these problematic observations
#' afterwards and therefore the final covariance matrix with imputed data is not
#' the same as the working covariance matrix (which may be based on prelminary mean
#' imputation).
#'
#' @param data a data frame or matrix with the data.
#' @param weights sampling weights.
#' @param outind an indicator vector for the outliers with \code{1} indicating
#' an outlier.
#' @param errors matrix of indicators for items which failed edits.
#' @param missing.matrix the missingness matrix can be given as input. Otherwise,
#' it will be recalculated.
#' @param alpha scalar giving the weight attributed to an item that is failing.
#' @param beta minimal overlap to accept a donor.
#' @param reweight.out if \code{TRUE}, the outliers are redefined.
#' @param c tuning constant when redefining the outliers (cutoff for Mahalanobis
#' distance).
#' @param preliminary.mean.imputation assume the problematic observation is at
#' the mean of good observations.
#' @param monitor if \code{TRUE} verbose output.
#' @return \code{POEM} returns a list whose first component \code{output} is a
#' sub-list with the following components:
#' \describe{
#'   \item{\code{preliminary.mean.imputation}}{Logical. \code{TRUE} if preliminary
#'   mean imputation should be used}
#'   \item{\code{completely.missing}}{Number of observations with no observed values}
#'   \item{\code{good.values}}{Weighted number of of good values (not missing, not
#'   outlying, not erroneous)}
#'   \item{\code{nonoutliers.before}}{Number of nonoutliers before reweighting}
#'   \item{\code{weighted.nonoutliers.before}}{Weighted number of nonoutliers
#'   before reweighting}
#'   \item{\code{nonoutliers.after}}{Number of nonoutliers after reweighting}
#'   \item{\code{weighted.nonoutliers.after}}{Weighted number of nonoutliers after
#'   reweighting}
#'   \item{\code{old.center}}{Coordinate means after weighting, before imputation}
#'   \item{\code{old.variances}}{Coordinate variances after weighting, before imputation}
#'   \item{\code{new.center}}{Coordinate means after weighting, after imputation}
#'   \item{\code{new.variances}}{Coordinate variances after weighting, after imputation}
#'   \item{\code{covariance}}{Covariance (of standardised observations) before imputation}
#'   \item{\code{imputed.observations}}{Indices of observations with imputated values}
#'   \item{\code{donors}}{Indices of donors for imputed observations}
#'   \item{\code{new.outind}}{Indices of new outliers}
#' }
#' The further component returned by \code{POEM} is:
#' \describe{
#'   \item{\code{imputed.data}}{Imputed data set}
#' }
#' @author Beat Hulliger
#' @references Béguin, C. and Hulliger B., (2002), EUREDIT Workpackage x.2
#' D4-5.2.1-2.C Develop and evaluate new methods for statistical outlier
#' detection and outlier robust multivariate imputation, Technical report,
#' EUREDIT 2002.
#' @examples
#' data(bushfirem, bushfire.weights)
#' outliers <- rep(0, nrow(bushfirem))
#' outliers[31:38] <- 1
#' imp.res <- POEM(bushfirem, bushfire.weights, outliers,
#' preliminary.mean.imputation = TRUE)
#' print(imp.res$output)
#' var(imp.res$imputed.data)
#' @export
#' @importFrom stats mahalanobis qchisq pnorm
POEM <- function(data, weights, outind, errors, missing.matrix, alpha = 0.5,
                 beta = 0.5, reweight.out = FALSE, c = 5,
                 preliminary.mean.imputation = FALSE, monitor = FALSE) {

  # ------- preparation -------

  # transform data to matrix
  if (!is.matrix(data)) {data <- as.matrix(data)}

  # number of rows
  n <- nrow(data)

  # number of columns
  p <- ncol(data)

  # set weights to 1 if missing
  if (missing(weights)) {weights <- rep(1, n)}

  # reformat 'outind' as numeric
  if (is.logical(outind)) {outind <- as.numeric(outind)}

  # if 'missing.matrix' is missing, create it
  if (missing(missing.matrix)) {missing.matrix <- (1 - is.na(data))}

  # if 'errors' is missing, create it
  if (missing(errors)) {errors <- matrix(1, nrow = n, ncol = p)}

  # initial tests (stops program if test fails)
  if (!is.vector(weights)) {stop("Weights not in vector form","\n")}
  if (!is.matrix(missing.matrix)) {stop("Missing values not in matrix form","\n")}
  if (!is.vector(outind)) {stop("outind not in vector form","\n")}
  if (!is.matrix(errors)) {stop("Errors not in matrix form","\n")}
  if (length(weights) != n) {stop("Wrong length of weights")}
  if (nrow(missing.matrix) != n | ncol(missing.matrix) != p) {
    stop("Missing values matrix do not have same dimensions as data","\n")
  }
  if (length(outind) != n) {stop("Wrong length of outind","\n")}
  if (nrow(errors) != n | ncol(errors) != p) {
    stop("Errors matrix do not have same dimensions as data","\n")
  }
  if (sum(is.na(weights)) > 0) {stop("Missing values in weights","\n")}
  if (sum(is.na(missing.matrix)) > 0) {stop("Missing values in missing data matrix","\n")}
  if (sum(is.na(outind)) > 0) {stop("Missing values in outind","\n")}
  if (sum(is.na(errors)) > 0) {stop("Missing values in errors matrix","\n")}

  # count number of completely missing observations
  comp.miss <- apply(is.na(data), 1, prod)
  cat("\n Number of completely missing observations ", sum(comp.miss), "\n")

  # start computation time
  calc.time <- proc.time()

  # set all missing values to zero
  old.data <- data
  data[!(missing.matrix)] <- 0

  # computation of alpha_ij
  alpha.ij <- missing.matrix * (alpha ^ (1 - errors))
  good.values <- apply(weights * alpha.ij, 2, sum)
  old.number.of.nonoutliers <- sum(1 - outind)
  old.weighted.sum.of.nonoutliers <- sum(weights * (1 - outind))
  missing.errors <- missing.matrix * errors

  # computation of center
  center <- apply((1 - outind) * weights * alpha.ij * data, 2, sum) /
    apply((1 - outind) * weights * alpha.ij, 2, sum)

  # centering of data
  data <- sweep(data, 2, center, "-")

  # computation of coordinates variances
  variances <- apply((1 - outind) * weights * alpha.ij * data ^ 2, 2, sum)
  variances <- variances / apply((1 - outind) * weights * alpha.ij, 2, sum)

  # exit function if there are zero variances
  if (sum(variances == 0) > 0) {
    zero.variances <- which(variances == 0)
	  cat("Warning: Variable(s)", zero.variances, "has (have) zero variance(s)\n")
	  stop("\nRemove these variables or reduce the set of outliers\n")
  }

  # standardization of data
  data <- sweep(data, 2, sqrt(variances), "/")

  # computation of covariance matrix
  covariance <- (t(alpha.ij * data) %*% ((1 - outind) * weights * alpha.ij * data))

  # compute covariance matrix
  if (!preliminary.mean.imputation) {

    # Formula (12) on page 121 of the Euredit Deliverable
	  covariance <- covariance / (t(alpha.ij) %*% ((1 - outind) * weights * alpha.ij))

	  # if preliminary.mean.imputation = FALSE, print warning if
	  # covariance matrix is positive definite
	  if (determinant(covariance)$sign < 0) {
	    cat("Warning: Covariance matrix not positive definite with original data including missing values\n")
	    cat("         Choose option preliminary.mean.imputation = TRUE\n")
	  }
  } else {

    # preliminary mean imputation is Formula (13) on page 122 of the Euredit Deliverable.
    # This is equivalent to setting tilde x to zero. And this in turn
    # is equivalent to setting x to the mean before standardising.
    covariance <- covariance / sum((1-outind)*weights)
  }

  # print progress to console
  if (monitor) {
    cat("Covariance matrix of standardised observations\n")
    print(covariance)
  }

  # reweighting of outliers
  if (reweight.out) {

    # compute Mahalanobis distances
    MD <- mahalanobis(alpha.ij * data, rep(0, p), covariance)
		MD <- p ^ 2 * MD / apply(alpha.ij, 1, sum) ^ 2

		# print warning if there are negative MD
		if (min(MD) < 0) {cat("Warning: Negative Mahalanobis distances\n")}

		# reweight outliers
		outind <- 1 - (MD > (c * qchisq(2 * pnorm(1) - 1, p)))

		# print progress to console
	  cat("New set of", sum(outind), "outliers generated\n")
  }

  # list of observations to be imputed
  observations.with.errors <- (apply(missing.errors, 1, sum) < p)
  to.be.imputed <- which(outind | observations.with.errors)
  imputed <- rep(0, n)

  # ------- Start of imputation process -------

  # list of complete and correct donors
  complete.donors <- apply((1 - outind) * missing.errors, 1, sum) == p

  # print progress to console
  cat("\n Number of complete and error-free observations: ", sum(complete.donors), "\n")

  # loop over observations to be imputed
  for (observation in to.be.imputed) {

    # indicator of potential donnors for the observation
    if (outind[observation] == 0 & beta < 1) {

      # conditions to be satisfied for indicator
      cond1 <- apply((1 - outind) * sweep(missing.errors, 2, (missing.errors[observation, ]), "*"), 1, sum) >= beta * p
      cond2 <- apply(sweep(1 - missing.errors, 2, (1 - missing.matrix[observation, ]), "*"), 1, sum) == 0
      cond3 <- apply(sweep(1 - missing.errors, 2, (1 - errors[observation, ]), "*"), 1, sum) == 0

      # find potential donors
      potential.donors <- cond1 & cond2 & cond3

    } else {

      potential.donors <- complete.donors

    }

    # exclude the observation as its own donor
    potential.donors[observation] <- FALSE

    # depending on the result of the first expression, first {} or second {} is evaluated
    switch(sum(potential.donors) + 1,
           {cat("\n No donor for observation ", observation, "\n")
             cat("All complete error-free observations used as donors.\n")
             potential.donors <- complete.donors},
           {cat("\n Only one donor for observation ", observation, "\n")})

    # print progress to console
    if (monitor) {cat("Number of potential donors ", sum(potential.donors), "\n")}

    # distances of the observation to potential donors
    distances.to.donors <-
      mahalanobis(sweep(data[potential.donors, ], 2, data[observation, ], "-") *
                    sweep(alpha.ij[potential.donors, ], 2,
                          alpha.ij[observation, ], "*"), rep(0, p), covariance)

    # transform distances
    distances.to.donors <-
      p ^ 2 * distances.to.donors /
      apply(sweep(alpha.ij[potential.donors, ], 2,
                  alpha.ij[observation, ], "*"), 1, sum) ^ 2

    # selection of the donor
    min.dist.to.donors <- min(distances.to.donors)

    # donors are the indices for a vector with the indices of the potential donors
    # potential.donors is an indicator vector of length n
    if (is.na(min.dist.to.donors)) {

      # if min. distance is missing
      donors <- 1:sum(potential.donors)

    } else {

      # if min. distance is negative
      if (min.dist.to.donors < 0) {

        # exit program
        stop("Warning: Minimal distance to nearest neighbour negative\n")

      } else {

        # get indices of donors with min. distance
        donors <- which(distances.to.donors == min.dist.to.donors)

      }
    }

    # one donor or multiple donors (in latter case sample)
    if (length(donors) == 1) {
      imputed[observation] <- which(potential.donors)[donors]
      } else {
        imputed[observation] <- which(potential.donors)[sample(donors, 1)]
      }

    # print progress to console
    if (monitor) {
      cat("Observation", observation, "imputed by donor ", imputed[observation], "\n")
      cat("distance to donor: ", min.dist.to.donors, "\n")
    }
  }

  # ------- IMPUTATION -------

  # create new data frame
  new.data <- old.data

  # imputation
  new.data[as.logical(outind), ] <- old.data[imputed[as.logical(outind)], ]
  non.outliers.errors <- which((1 - outind) & observations.with.errors)
  new.data[non.outliers.errors, ][missing.errors[non.outliers.errors,] < 1]  <-
    old.data[imputed[non.outliers.errors], ][missing.errors[non.outliers.errors, ] < 1]

  # compute new center and variances (with imputed data)
  new.center <- apply(weights * new.data, 2, sum) / sum(weights)
  new.variances <-
    apply(weights * sweep(new.data, 2, new.center, "-") ^ 2, 2, sum) / sum(weights)

	# ------- results -------

  # stop computation time
  calc.time <- proc.time() - calc.time

	# prepare output
  POEM.r <- list(
    preliminary.mean.imputation = preliminary.mean.imputation,
    completely.missing = sum(comp.miss),
    good.values = good.values,
		nonoutliers.before = old.number.of.nonoutliers,
		weighted.nonoutliers.before = old.weighted.sum.of.nonoutliers,
		number.of.nonoutliers.after.reweighting = sum(1 - outind),
		weighted.number.of.nonoutliers.after.reweighting = sum(weights * (1 - outind)),
		old.center = center,
		old.variances = variances,
    new.center = new.center,
		new.variances = new.variances,
    covariance = covariance,
		imputed.observations = to.be.imputed,
		donors = imputed[to.be.imputed],
		outind = outind)

	# return output
	return(invisible(list(output = POEM.r, imputed.data = new.data)))

	# output to console
	cat("\n", "POEM has imputed", length(to.be.imputed),
	    "observations(s) in", calc.time, "seconds.", "\n", "\n")
}

