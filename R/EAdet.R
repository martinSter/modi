#' Epidemic Algorithm for detection of multivariate outliers in incomplete survey data
#'
#' In \code{EAdet} an epidemic is started at a center of the data. The epidemic
#' spreads out and infects neighbouring points (probabilistically or deterministically).
#' The last points infected are outliers. After running \code{EAdet} an imputation
#' with \code{EAimp} may be run.
#'
#' The form and parameters of the transmission function should be chosen such that the
#' infection times have at least a range of 10. The default cutting point to decide on
#' outliers is the median infection time plus three times the mad of infection times.
#' A better cutpoint may be chosen by visual inspection of the cdf of infection times.
#' \code{EAdet} calls the function \code{EA.dist}, which passes the counterprobabilities
#' of infection (a \eqn{n * (n - 1) / 2} size vector!) and three parameters (sample
#' spatial median index, maximal distance to nearest neighbor and transmission distance =
#' reach) as arguments to \code{EAdet}. The distances vector may be too large to be passed
#' as arguments. Then either the memory size must be increased. Former versions of the
#' code used a global variable to store the distances in order to save memory.
#'
#' @param data a data frame or matrix with data.
#' @param weights a vector of positive sampling weights.
#' @param reach if \code{reach = "max"} the maximal nearest neighbor distance is
#' used as the basis for the transmission function, otherwise the weighted
#' \eqn{(1 - (p + 1) / n)} quantile of the nearest neighbor distances is used.
#' @param transmission.function form of the transmission function of distance d:
#' \code{"step"} is a heaviside function which jumps to \code{1} at \code{d0},
#' \code{"linear"} is linear between \code{0} and \code{d0}, \code{"power"} is
#' \code{(beta*d+1)^(-p)} for \code{p = ncol(data)} as default, \code{"root"} is the
#' function \code{1-(1-d/d0)^(1/maxl)}.
#' @param power sets \code{p = power}.
#' @param distance.type distance type in function \code{dist()}.
#' @param maxl maximum number of steps without infection.
#' @param plotting if \code{TRUE}, the cdf of infection times is plotted.
#' @param monitor if \code{TRUE}, verbose output on epidemic.
#' @param prob.quantile if mads fail, take this quantile absolute deviation.
#' @param random.start if \code{TRUE}, take a starting point at random instead of the
#' spatial median.
#' @param fix.start force epidemic to start at a specific observation.
#' @param threshold infect all remaining points with infection probability above
#' the threshold \code{1-0.5^(1/maxl)}.
#' @param deterministic if \code{TRUE}, the number of infections is the expected
#' number and the infected observations are the ones with largest infection probabilities.
#' @param rm.missobs set \code{rm.missobs=TRUE} if completely missing observations
#' should be discarded. This has to be done actively as a safeguard to avoid mismatches
#' when imputing.
#' @param verbose more output with \code{verbose=TRUE}.
#' @return \code{EAdet} returns a list whose first component \code{output} is a sub-list
#' with the following components:
#' \describe{
#'   \item{\code{sample.size}}{Number of observations}
#'   \item{\code{discarded.observations}}{Indices of discarded observations}
#'   \item{\code{missing.observations}}{Indices of completely missing observations}
#'   \item{\code{number.of.variables}}{Number of variables}
#'   \item{\code{n.complete.records}}{Number of records without missing values}
#'   \item{\code{n.usable.records}}{Number of records with less than half of values
#'   missing (unusable observations are discarded)}
#'   \item{\code{medians}}{Component wise medians}
#'   \item{\code{mads}}{Component wise mads}
#'   \item{\code{prob.quantile}}{Use this quantile if mads fail, i.e. if one of the mads is 0}
#'   \item{\code{quantile.deviations}}{Quantile of absolute deviations}
#'   \item{\code{start}}{Starting observation}
#'   \item{\code{transmission.function}}{Input parameter}
#'   \item{\code{power}}{Input parameter}
#'   \item{\code{maxl}}{Maximum number of steps without infection}
#'   \item{\code{min.nn.dist}}{Maximal nearest neighbor distance}
#'   \item{\code{transmission.distance}}{\code{d0}}
#'   \item{\code{threshold}}{Input parameter}
#'   \item{\code{distance.type}}{Input parameter}
#'   \item{\code{deterministic}}{Input parameter}
#'   \item{\code{number.infected}}{Number of infected observations}
#'   \item{\code{cutpoint}}{Cutpoint of infection times for outlier definition}
#'   \item{\code{number.outliers}}{Number of outliers}
#'   \item{\code{outliers}}{Indices of outliers}
#'   \item{\code{duration}}{Duration of epidemic}
#'   \item{\code{computation.time}}{Elapsed computation time}
#'   \item{\code{initialisation.computation.time}}{Elapsed computation time for
#'   standardisation and calculation of distance matrix}
#' }
#' The further components returned by \code{EAdet} are:
#' \describe{
#'   \item{\code{infected}}{Indicator of infection}
#'   \item{\code{infection.time}}{Time of infection}
#'   \item{\code{outind}}{Indicator of outliers}
#' }
#' @author Beat Hulliger
#' @references Béguin, C. and Hulliger, B. (2004) Multivariate outlier detection in
#' incomplete survey data: the epidemic algorithm and transformed rank correlations,
#' JRSS-A, 167, Part 2, pp. 275-294.
#' @seealso \code{\link{EAimp}} for imputation with the Epidemic Algorithm.
#' @examples
#' data(bushfirem, bushfire.weights)
#' det.res <- EAdet(bushfirem, bushfire.weights)
#' print(det.res$output)
#' @export
#' @importFrom stats rbinom
#' @importFrom graphics plot abline
#' @importFrom utils memory.size
EAdet <- function(data, weights, reach = "max", transmission.function = "root",
                  power = ncol(data), distance.type = "euclidean", maxl = 5,
                  plotting = TRUE, monitor = FALSE, prob.quantile = 0.9,
                  random.start = FALSE, fix.start, threshold = FALSE,
                  deterministic = TRUE, rm.missobs = FALSE, verbose = FALSE) {

  # ------- preparation -------

  # transform data to matrix
  if (!is.matrix(data)) {data <- as.matrix(data)}

  # number of rows
  n <- nrow(data)

  # number of columns
  p <- ncol(data)

  # set weights to 1 if missing
  if (missing(weights)) {weights <- rep(1, n)}

  # finding the unit(s) with all items missing
  new.indices <- which(apply(is.na(data), 1, prod) == 0)
  miss.indices <- apply(is.na(data), 1, prod) == 1
  missobs <- which(miss.indices)
  discarded <- NA
  nfull <- n

  # removing the unit(s) with all items missing
  if ((length(new.indices) < n) & rm.missobs) {
    discarded <- missobs
    cat("Warning: missing observations", discarded, "removed from the data\n")
    data <- data[new.indices, ]
    weights <- weights[new.indices]
    n <- nrow(data)
  }

  # find complete and usable (valid information for at least half of var.) records
	complete.records <- apply(!is.na(data), 1, prod)
	usable.records <- apply(!is.na(data), 1, sum) >= p / 2

	# print progress to console
	if (verbose) {
	  cat("\n Dimensions (n,p):", n, p)
	  cat("\n Number of complete records ", sum(complete.records))
	  cat("\n Number of records with maximum p/2 variables missing ",
	      sum(usable.records), "\n")
	}

	# transform parameter power to double
  power <- as.single(power)

  # standardization of weights to sum to sample size
	np <- sum(weights)
	weights <- as.single((n * weights) / np)

	# start computation time
	calc.time <- proc.time()[1]

	# ------- calibraton and setup -------

	# compute medians
	medians <- apply(data, 2, weighted.quantile, w = weights, prob = 0.5)

	# sweep median from data
	data <- sweep(data, 2, medians, "-")

	# compute median absolute deviations
	mads <- apply(abs(data), 2, weighted.quantile, w = weights, prob = 0.5)

	# compute quantile absolute deviations
	qads <- apply(abs(data), 2, weighted.quantile, w = weights, prob = prob.quantile)

	# standardization
	if(sum(mads == 0) > 0) {

	  # output to console if some mads are 0
		cat("\n Some mads are 0. Standardizing with ", prob.quantile,
		    " quantile absolute deviations!")

		if(sum(qads == 0) > 0) {

		  # if some qads are 0, no standardization
		  cat("\n Some quantile absolute deviations are 0. No standardization!")

		} else {

		  # if no qads are 0, standardize with qads
		  data <- sweep(data, 2, qads, "/")

		}

	} else {

	  # if no mads are 0, standardize with mads
	  data <- sweep(data, 2, mads, "/")

	}

	# calculation of distances
	EA.dist.res <- EA.dist(data, n = n, p = p, weights = weights, reach = reach,
	                        transmission.function = transmission.function,
	                        power = power, distance.type = distance.type,
	                        maxl = maxl)

	# print progress to console
	if (monitor) {cat("\n\n Distances finished")}

	# print progress to console
	# The distances calculated by EA.dist are the
	# counterprobabilities in single precision.
  if (monitor) {
    cat("\n memory use after distances: ", memory.size())
    cat("\n Index of sample spatial median is ", EA.dist.res$output[1])
    cat("\n Maximal distance to nearest neighbor is ", EA.dist.res$output[2])
    cat("\n Transmission distance is ", EA.dist.res$output[3], "\n")
  }

	# ------- initialisation -------

	# print progress to console
	if (verbose) {cat("\n\n Initialisation of epidemic")}

	# initialization time
	comp.time.init <- proc.time()[1] - calc.time

	# print initialization time to console
	if(monitor) {cat("\n Initialisation time is ", comp.time.init)}

	# define starting point of infection
	if(random.start) {

	  # random starting point
	  start.point <- sample(1:n, 1, prob = weights)

	} else {

	  if(!missing(fix.start)) {

	    # if fix.start is defined, then this is the starting point
	    start.point <- fix.start

	  } else {

	    # else start with first obs.
	    start.point <- EA.dist.res$output[1]

	  }
	}

	# set time to 1
	time <- 1

	# initialize infected vector with FALSE
	infected <- rep(FALSE, n)

	# set starting point of inection to TRUE
	infected[c(start.point)] <- TRUE

  # initialize various things
	new.infected <- infected
	n.infected <- sum(infected)
	hprod <- rep(1, n)

	infection.time <- rep(0, n)
	infection.time[c(start.point)] <- time

	# ------- main loop -------
	repeat {

	  # print progress to console
		if (monitor) {cat("\n time = ", time, " , infected = ", n.infected)}

		time <- time + 1
		old.infected <- infected

		if(sum(new.infected) > 1) {

		  hprod[!infected] <-
		    hprod[!infected] * apply(sweep(
		      sweep(matrix(EA.dist.res$counterprobs[apply(
		        as.matrix(which(!infected)), 1, ind.dijs,
		        js = which(new.infected), n = n)], sum(new.infected),
		        n - n.infected), 1, weights[new.infected], "^"), 2,
		      weights[!infected], "^"), 2, prod)

		} else {

			if(sum(new.infected) == 1) {

			  hprod[!infected] <-
			    hprod[!infected] * EA.dist.res$counterprobs[apply(
			      as.matrix(which(!infected)), 1, ind.dijs,
			      js = which(new.infected), n = n)] ^ (weights[new.infected] *
			                                             weights[!infected])

			}

		}

		if (deterministic) {

			n.to.infect <- sum(1 - hprod[!infected]) # HRK: expected number of infections

			# do maxl trials for very small inf. prob.
			if (n.to.infect < 0.5) {

			  n.to.infect <- sum(1 - hprod[!infected] ^ maxl)

			}

			n.to.infect <- round(n.to.infect)

			infected[!infected] <-
			  rank(1 - hprod[!infected]) >= n - n.infected - n.to.infect

		} else {

			if (threshold) {

			  infected[!infected] <- hprod[!infected] <= 0.5 ^ (1 / maxl)

			} else {

			  infected[!infected] <-
			    as.logical(rbinom(n - n.infected, 1, 1 - hprod[!infected]))

			}

		}

		new.infected <- infected & (!old.infected)
		n.infected <- sum(infected)
		infection.time[new.infected] <- time

		# if all are infected, stop loop
		if(n.infected == n) {break}

		# if max. infection steps is reached, stop loop
		if((time - max(infection.time)) > maxl) {break}

		# start next iteration of loop
		next
	}

	# duration of infection
	duration <- max(infection.time)

	# print progress to console
	if (verbose) {cat("\n memory use after epidemic: ", memory.size())}

	# stop computation time
	calc.time <- round(proc.time()[1] - calc.time, 5)

  # default cutpoint
  med.infection.time <- weighted.quantile(infection.time, weights, 0.5)
	mad.infection.time <-
	  weighted.quantile(abs(infection.time - med.infection.time), weights, 0.5)

	# print progress to console
	if (verbose) {
	  cat("\n med and mad of infection times: ", med.infection.time,
	      " and ", mad.infection.time)
	}

	if (mad.infection.time == 0) {
	  mad.infection.time <- med.infection.time
	}

	cutpoint <- min(med.infection.time + 3 * mad.infection.time, duration)

	# print progress to console
	if (verbose) {cat("\n Proposed cutpoint is ", min(cutpoint, duration))}

  # blowing up to full length
  infectedn <- logical(n)
  infectedn[infected] <- TRUE
  infectednfull <- logical(nfull)

  if (nfull > n) {

    infectednfull[new.indices] <- infectedn

    }	else {

      infectednfull <- infectedn

    }

  # initialize empty vector for infection times
  inf.time <- rep(NA, nfull)

  # get infection times
  if (nfull > n) {

    inf.time[new.indices] <- infection.time

    } else {

      inf.time <- infection.time

    }

  # set infection time of completely missing to NA
  inf.time[!infectednfull] <- NA

  # outliers full sample
  outlier <- (inf.time >= cutpoint)

  # get indices of outliers
  outlier.ind <- which(outlier)

  # ------- plotting -------
  # not infected are set to high inf.time to show better
  # the healthy on a graph of infection times

  if(plotting) {
    plot.time <- inf.time
    plot.time[!infectednfull] <- ceiling(1.2 * duration)
    ord <- order(plot.time)
    plot(plot.time[ord], cumsum(weights[ord]), xlab = "infection time",
         ylab = "(weighted) cdf of infection time")
    abline(v = cutpoint)
  }

  # ------- results -------

  # prepare output
	EAdet.r <- list(
	  sample.size = n, discarded.observations = discarded,
	  missing.observations = missobs, number.of.variables = p,
	  n.complete.records = sum(complete.records),
	  n.usable.records = sum(usable.records), medians = medians,
	  mads = mads, prob.quantile = prob.quantile,
	  quantile.deviations = qads, start = start.point,
	  transmission.function = transmission.function, power = power,
	  maxl = maxl, max.min.di = EA.dist.res$output[2],
	  transmission.distance = EA.dist.res$output[3], threshold = threshold,
	  distance.type = distance.type, deterministic = deterministic,
	  number.infected = n.infected, cutpoint = cutpoint,
	  number.outliers = sum(outlier), outliers = outlier.ind,
	  duration = duration, computation.time = calc.time,
	  initialisation.computation.time = comp.time.init)

  # output to console
	cat("\n", "EA detection has finished with", n.infected, "infected points in",
	    calc.time[1], "seconds.")

  # return output
  return(invisible(list(output = EAdet.r, infected = infectednfull,
                        infection.time = inf.time, outind = outlier)))
}

