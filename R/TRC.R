#' Transformed rank correlations for multivariate outlier detection
#'
#' \code{TRC} starts from bivariate Spearman correlations and obtains
#' a positive definite covariance matrix by back-transforming robust
#' univariate medians and mads of the eigenspace. \code{TRC} can cope
#' with missing values by a regression imputation using the a robust
#' regression on the best predictor and it takes sampling weights
#' into account.
#'
#' \code{TRC} is similar to a one-step OGK estimator where the starting
#' covariances are obtained from rank correlations and an ad hoc missing
#' value imputation plus weighting is provided.
#'
#' @param data a data frame or matrix with the data.
#' @param weights sampling weights.
#' @param overlap minimum number of jointly observed values for calculating
#' the rank correlation.
#' @param mincor minimal absolute correlation to impute.
#' @param robust.regression type of regression: \code{"irls"} is iteratively
#' reweighted least squares M-estimator, \code{"rank"} is based on the rank
#' correlations.
#' @param gamma minimal number of jointly observed values to impute.
#' @param prob.quantile if mads are 0, try this quantile of absolute deviations.
#' @param alpha \code{(1 - alpha)} Quantile of F-distribution is used for cut-off.
#' @param md.type type of Mahalanobis distance when missing values occur:
#' \code{"m"} marginal (default), \code{"c"} conditional.
#' @param monitor if \code{TRUE}, verbose output.
#' @return \code{TRC} returns a list whose first component \code{output} is a
#' sublist with the following components:
#' \describe{
#'   \item{\code{sample.size}}{Number of observations}
#'   \item{\code{number.of.variables}}{Number of variables}
#'   \item{\code{number.of.missing.items}}{Number of missing values}
#'   \item{\code{significance.level}}{\code{1 - alpha}}
#'   \item{\code{computation.time}}{Elapsed computation time}
#'   \item{\code{medians}}{Componentwise medians}
#'   \item{\code{mads}}{Componentwise mads}
#'   \item{\code{center}}{Location estimate}
#'   \item{\code{scatter}}{Covariance estimate}
#'   \item{\code{robust.regression}}{Input parameter}
#'   \item{\code{md.type}}{Input parameter}
#'   \item{\code{cutpoint}}{The default threshold MD-value for the cut-off of outliers}
#' }
#' The further components returned by \code{TRC} are:
#' \describe{
#'   \item{\code{outind}}{Indicator of outliers}
#'   \item{\code{dist}}{Mahalanobis distances (with missing values)}
#' }
#' @author Beat Hulliger
#' @references Béguin, C. and Hulliger, B. (2004) Multivariate outlier detection in
#' incomplete survey data: the epidemic algorithm and transformed rank correlations,
#' JRSS-A, 167, Part 2, pp. 275-294.
#' @examples
#' data(bushfirem, bushfire.weights)
#' det.res <- TRC(bushfirem, weights = bushfire.weights)
#' PlotMD(det.res$dist, ncol(bushfirem))
#' print(det.res)
#' @export
#' @importFrom stats qnorm mahalanobis qf median
TRC <- function(data, weights, overlap = 3, mincor = 0,
                robust.regression = "rank", gamma = 0.5,
                prob.quantile = 0.75, alpha = 0.05,
                md.type = "m", monitor = FALSE) {

  # start computation time
  calc.time <- proc.time()

  # ------- preprocessing -------

  # transform data to matrix
	if(!is.matrix(data)) {data <- as.matrix(data)}

  # number of rows
	n <- nrow(data)

	# number of columns
	p <- ncol(data)

	# correct alpha if alpha smaller 0.5
	if (alpha < 0.5) {alpha <- 1 - alpha}

	# set weights to 1 if missing
	if (missing(weights)) {(weights <- rep(1, n))}

  # stop program if weights do have wrong dimension
	if (length(weights) != n) {stop("\n Sampling weights of wrong length")}

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

	# create matrix with missings
	missing.matrix <- 1 - is.na(data)

	# output number of missing variables
	cat("\n Number of missing items: ", sum(is.na(data)),
	    ", percentage of missing items: ", mean(is.na(data)), " \n")

	# print progress to console
	if (monitor) {
	  cat("End of preprocessing in", proc.time() - calc.time, "seconds \n")
	}

	# ------- estimation of location and scatter -------

  # starting time of estimation of location and scatter
	if (monitor) {spearman.time <- proc.time()}

  # compute weighted medians
	medians <- apply(data, 2, weighted.quantile, w = weights)

	# correction constant for median absolute deviations
	correction.constant.mad <- 1 / qnorm(0.75)

	# compute median absolute deviations
	mads <- correction.constant.mad *
	  apply(abs(sweep(data, 2, medians, "-")), 2, weighted.quantile, w = weights)

  # deal with observations with 'mads' smaller or equal to 0
	if(sum(mads <= 0) > 0) {

	  # use 'prob.quantile'
		cat("Some mads are 0. Using", prob.quantile, "quantile absolute deviations!\n")
		mads <- (1 / qnorm(0.5 * (1 + prob.quantile))) *
		  apply(abs(sweep(data, 2, medians, "-")), 2, weighted.quantile,
		        w = weights, prob = prob.quantile)

		# stop program if there are still mads' smaller or equal to 0
		if(sum(mads <= 0) > 0) {
			cat("The following variable(s) have", prob.quantile,
			    "quantile absolute deviations equal to 0 :", which(mads == 0), "\n")
			stop("Remove these variables or increase the quantile probablity\n")
		}
	}

	# if no missing values
	if (prod(missing.matrix) == 1) {

	  # initialize matrix
		weighted.ranks <- matrix(0, n, p)

		# loop over columns
		for (i in 1:p) {

			weighted.ranks[ ,i] <-
			  (apply(data[ , i, drop = FALSE], 1, weightsum,
			         weights = weights, observations = data[ ,i]) +
			     0.5 * apply(data[ , i, drop = FALSE], 1, weightsum,
			                 weights = weights, observations = data[ ,i], lt = FALSE) +
			     0.5)
		}

		# compute scatter
		scatter <- (12 * (t(weighted.ranks) %*% (weights * weighted.ranks)) /
		              (t(missing.matrix) %*% (weights * missing.matrix))^3 - 3)

		# set scatter to 1 where it is larger than 1
		scatter[scatter > 1] <- 1

		# set scatter to -1 where it is smaller than -1
		scatter[scatter < (-1)] <- -1

		# standardization
		scatter <- 2 * sin(pi * scatter / 6)

    # print progress to console
		if (monitor) {
			cat("Spearman Rank Correlations (truncated and standardized):\n")
			print(scatter)
			cat("End of Spearman rank correlations estimations in",
			    proc.time() - spearman.time, "seconds\n")
		}

    # print progress to console
		if (monitor) {cat("No imputation\n")}

		# if missing values (very long else statement)
	} else {

	  # initialize matrix
		scatter <- matrix(0, p, p)

		# ???
		size.of.cor.sets <- t(missing.matrix) %*% missing.matrix

		# print warning
		if (sum(size.of.cor.sets < overlap) > 0) {
			cat("Warning: ", (sum(size.of.cor.sets < overlap) - p) / 2,
			    " couples of variables have less than ", overlap,
			    " observations in common, therefore their rank
			    correlations will be set to 0.\n")
		}

		# print progress to console
		if (monitor) {cat("Computing Spearman Rank Correlations :\n")}

    # start loop
		for (i in 1:(p - 1)) {

		  # print progress to console
			if (monitor) {cat("i=", i, "\n")}

		  # start inner loop
			for (j in (i + 1):p) {

			  # print progress to console
				if (monitor) {cat(" j=", j, "\n")}

				if (size.of.cor.sets[i,j] >= overlap) {

					common.observations <- missing.matrix[ ,i] & missing.matrix[ ,j]

					weighted.ranks.i <- (apply(data[common.observations, i, drop = FALSE], 1,
					                           weightsum, weights = weights[common.observations],
					                           observations = data[common.observations,i]) +
					                       0.5 * apply(data[common.observations, i, drop = FALSE], 1,
					                                   weightsum,
					                                   weights = weights[common.observations],
					                                   observations = data[common.observations,i],
					                                   lt = FALSE) + 0.5)

					weighted.ranks.j <- (apply(data[common.observations, j, drop = FALSE], 1,
					                           weightsum, weights = weights[common.observations],
					                           observations = data[common.observations,j]) +
					                       0.5 * apply(data[common.observations, j, drop = FALSE], 1,
					                                   weightsum,
					                                   weights = weights[common.observations],
					                                   observations = data[common.observations,j],
					                                   lt = FALSE) + 0.5)

					scatter[i,j] <-
					  12 * sum(weights[common.observations] * weighted.ranks.i * weighted.ranks.j) /
					  sum(weights[common.observations])^3 - 3
				}
			}
		}

		# compute scatter
		scatter <- scatter + t(scatter) + diag(p)

		# set scatter to 1 where it is larger than 1
		scatter[scatter > 1] <- 1

		# set scatter to -1 where it is smaller than -1
		scatter[scatter < (-1)] <- -1

		# standardization
		scatter <- 2 * sin(pi * scatter / 6)

		# print progress to console
		if (monitor) {
    cat("Spearman Rank Correlations (truncated and standardized):\n")
		print(scatter)
		cat("End of Spearman rank correlations estimations in",
		    proc.time() - spearman.time, "seconds\n")
		}

		# ------- ad hoc imputation of missing values -------

		# start imputation time
		imputation.time <- proc.time()

		# get variables to be imputed
		variables.to.be.imputed <- which(apply(missing.matrix, 2, prod) == 0)

		# regressors correlations with small support are set to 0
    regressors.cor <- (scatter - diag(p))[variables.to.be.imputed, ] *
      (size.of.cor.sets[variables.to.be.imputed, ] >= (gamma * n))

		regressors.cor <- as.matrix(t(regressors.cor))

		# print progress to console
		if (monitor) {cat("Regressors correlations\n", regressors.cor)}

		# if there are more than one variables to be imputed
		if (length(variables.to.be.imputed) > 1) {

		  regressors.list.ordered <- apply(-abs(regressors.cor), 1, order)

		} else {

		  regressors.list.ordered <- as.matrix(t(order(-abs(regressors.cor))))

		}

    # loop over variables to be imputed
		for (v in 1:length(variables.to.be.imputed)) {

		  observations.to.be.imputed  <- (!missing.matrix[ ,variables.to.be.imputed[v]])

		  # print progress to console
		  if (monitor) {cat("Variable", variables.to.be.imputed[v], ":\n")}

      # loop over candidate regressors
      r <- 0

			repeat {

			  r <- r + 1

        # output warning and stop program
			  if (abs(regressors.cor[v, regressors.list.ordered[r,v]]) < mincor) {

			    # first in regressors.list.ordered is largest cor (output warning)
					cat("No eligible regressor found for variable", v,
					    "observation(s)", which(observations.to.be.imputed),
					    ".\n Try to relax the regressor eligibility conditions.\n")

					stop()
			  }

			  # if condition is satisfied, skip current iteration
				if (sum(observations.to.be.imputed &
				        missing.matrix[ ,regressors.list.ordered[r,v]]) == 0) {next}


				observations.imputed.by.r.on.v <-
				  which(observations.to.be.imputed &
				          missing.matrix[,regressors.list.ordered[r,v]])


				k <- length(observations.imputed.by.r.on.v)


				common.observations <- missing.matrix[ ,variables.to.be.imputed[v]] &
				  missing.matrix[ ,regressors.list.ordered[r,v]]

				# robust regression depending on type of regression
				if (robust.regression == "irls") {

					regression.coeff <-
					  MASS::rlm(data[common.observations, variables.to.be.imputed[v]] ~
					              data[common.observations, regressors.list.ordered[r,v]],
					            weights = weights[common.observations])$coefficients

				} else {

					regression.coeff <-
					  c(0, regressors.cor[v, regressors.list.ordered[r,v]] *
					      mads[variables.to.be.imputed[v]] / mads[regressors.list.ordered[r,v]])

					regression.coeff[1] <- medians[variables.to.be.imputed[v]] -
					  regression.coeff[2] * medians[regressors.list.ordered[r,v]]

				}


				data[observations.imputed.by.r.on.v, variables.to.be.imputed[v]] <-
				  matrix(c(rep(1, k), data[observations.imputed.by.r.on.v,
				                           regressors.list.ordered[r,v]]), k, 2) %*%
				  regression.coeff

				observations.to.be.imputed[observations.imputed.by.r.on.v] <- FALSE

				# print progress to console
				if (monitor) {
				  cat(" ", k, "observations imputed using regressor",
				      regressors.list.ordered[r,v], "(cor=",
				      scatter[variables.to.be.imputed[v], regressors.list.ordered[r,v]],
				      "slope =", regression.coeff[2], "intercept =",
				      regression.coeff[1], ")\n")
				}

				# if condition is satisfied, skip current iteration
				if (sum(observations.to.be.imputed) > 0) {next}

				# stop iteration
				break

			} # end of repeat loop

		} # end of for loop

		# print progress to console
		if (monitor) {cat("End of imputation in", proc.time()-imputation.time, "seconds\n")}

	} # end of else statement

	# scatter
	scatter <- t(t(mads * scatter) * mads)
	new.basis <- eigen(scatter)$vectors
	data <- data %*% new.basis

	# center
	center <- apply(data, 2, weighted.quantile, w = weights)

	# scatter
	scatter <-
	  (correction.constant.mad *
	     apply(abs(sweep(data, 2, center, "-")), 2, weighted.quantile, w = weights))^2

	# if some scatter elements are 0
	if(sum(scatter == 0) > 0) {

	  # show warning
		cat("Some mads are 0. Using", prob.quantile, "quantile absolute deviations!\n")

	  # scatter
	  scatter <-
	    ((1 / qnorm(0.5 * (1 + prob.quantile))) * apply(abs(sweep(data, 2, center, "-")),
	                                                    2, weighted.quantile, w = weights,
	                                                    prob = prob.quantile))^2

	  # stop program if there are still 0 scatter elements
		if(sum(scatter == 0) > 0) {
			stop("Please, increase the quantile probability\n")
		}
	}

	# get center and scatter of data
	center <- as.vector(new.basis %*% center)
	scatter <- as.matrix(new.basis %*% diag(scatter) %*% t(new.basis))
	data <- data %*% t(new.basis)

	# ------- Mahalanobis distances -------

  # set data to missing where missing
	data[!missing.matrix] <- NA

	# get patterns
	s.patterns <-
	  apply(matrix(as.integer(is.na(data)), n, p), 1, paste, sep = "", collapse = "")

	# order missingness patterns
	perm <- order(s.patterns)

	# order data and patterns according to 'perm'
	data <- data[perm, ]
	s.patterns <- s.patterns[perm]

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
	indices <- (!missing.items[1, ])

  # type of Mahalanobis distance
	if (md.type == "c") {
	   metric <- solve(scatter)
	  } else {
	    metric <- scatter
	  }

	# compute distances
	dist <- mahalanobis(data[1:s.id[1], indices, drop = FALSE], center[indices],
	                    metric[indices,indices],
	                    inverted = (md.type == "c")) * p / (p - nb.missing.items[1])

	# ???
	if (S > 1) {

	  # start loop
	  for (i in 2:S) {
	    indices <- (!missing.items[i, ])
	    dist <-
	      c(dist, mahalanobis(data[(s.id[i - 1] + 1):s.id[i], indices, drop = FALSE],
	                          center[indices], metric[indices, indices, drop = FALSE],
	                          inverted = (md.type == "c")) * p / (p - nb.missing.items[i]))
	  }
	}

	# ------- choice of outliers -------

	# nominate the outliers using the original numbering (without discarded obs.)
	cutpoint <- qf(alpha, p, n - p) / qf(0.5, p, n - p) * median(dist)
	dist <- dist[order(perm)]
	good <- (1:n)[dist < cutpoint]
	outliers <- (1:n)[-good]

	# outliers and distances in original numbering with full dataset
  outn <- logical(n)
  outn[outliers] <- TRUE
  outnfull <- logical(nfull)
  outnfull[new.indices] <- outn
  distnfull <- rep(NA, nfull)
  distnfull[new.indices] <- dist

	# ------- results -------

	# stop computation time
	calc.time <- proc.time() - calc.time

	# prepare output
	res <- list(
	  sample.size = n,
	  number.of.variables = p,
	  number.of.missing.items = nb.missing.items,
	  significance.level = alpha,
	  computation.time = calc.time,
	  medians = medians,
	  mads = mads,
	  center = center,
	  scatter = scatter,
	  robust.regression = robust.regression,
	  md.type = md.type,
	  cutpoint = cutpoint,
	  outind = outnfull,
	  dist = distnfull
	  )

	# output to console
	message(paste0("\n", "TRC has detected ", length(outliers), " outlier(s) in ",
	        round(calc.time[1], 2), " seconds.\n"))

	class(res) <- "TRC.r"
	res
}

