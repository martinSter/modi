#' Epidemic Algorithm for imputation of multivariate outliers in incomplete
#' survey data.
#'
#' After running \code{EAdet} an imputation of the detected outliers with
#' \code{EAimp} may be run.
#'
#' \code{EAimp} uses the distances calculated in \code{EAdet} (actually the
#' counterprobabilities, which are stored in a global data set) and starts an
#' epidemic at each observation to be imputed until donors for the missing values
#' are infected. Then a donor is selected randomly.
#'
#' @param data a data frame or matrix with the data.
#' @param weights a vector of positive sampling weights.
#' @param outind a logical vector with component \code{TRUE} for outliers.
#' @param reach reach of the threshold function (usually set to the maximum
#' distance to a nearest neighbour, see internal function \code{.EA.dist}).
#' @param transmission.function form of the transmission function of distance d:
#' \code{"step"} is a heaviside function which jumps to \code{1} at \code{d0},
#' \code{"linear"} is linear between \code{0} and \code{d0}, \code{"power"} is
#' \code{beta*d+1^(-p)} for \code{p=ncol(data)} as default, \code{"root"} is the
#' function \code{1-(1-d/d0)^(1/maxl)}.
#' @param power sets \code{p=power}, where \code{p} is the parameter in the above
#' transmission function.
#' @param distance.type distance type in function \code{dist()}.
#' @param maxl maximum number of steps without infection.
#' @param monitor if \code{TRUE} verbose output on epidemic.
#' @param threshold Infect all remaining points with infection probability above
#' the threshold \code{1-0.5^(1/maxl)}.
#' @param deterministic if \code{TRUE} the number of infections is the expected
#' number and the infected observations are the ones with largest infection
#' probabilities.
#' @param duration the duration of the detection epidemic.
#' @param kdon the number of donors that should be infected before imputation.
#' @param fixedprop if \code{TRUE} a fixed proportion of observations is infected
#' at each step.
#' @return \code{EAimp} returns a list with two components: \code{parameters} and
#' \code{imputed.data}.
#' \code{parameters} contains the following elements:
#' \describe{
#'   \item{\code{sample.size}}{Number of observations}
#'   \item{\code{number.of.variables}}{Number of variables}
#'   \item{\code{n.complete.records}}{Number of records without missing values}
#'   \item{\code{n.usable.records}}{Number of records with less than half of values
#'   missing (unusable observations are discarded)}
#'   \item{\code{duration}}{Duration of epidemic}
#'   \item{\code{reach}}{Transmission distance (\code{d0})}
#'   \item{\code{threshold}}{Input parameter}
#'   \item{\code{deterministic}}{Input parameter}
#'   \item{\code{computation.time}}{Elapsed computation time}
#' }
#' \code{imputed.data} contains the imputed data.
#' @author Beat Hulliger
#' @references Béguin, C. and Hulliger, B. (2004) Multivariate outlier detection in
#' incomplete survey data: the epidemic algorithm and transformed rank correlations,
#' JRSS-A, 167, Part 2, pp. 275-294.
#' @seealso \code{\link{EAdet}} for outlier detection with the Epidemic Algorithm.
#' @examples
#' data(bushfirem, bushfire.weights)
#' det.res <- EAdet(bushfirem, bushfire.weights)
#' imp.res <- EAimp(bushfirem, bushfire.weights, outind = det.res$outind,
#' reach = det.res$output$max.min.di, kdon = 3)
#' print(imp.res$output)
#' @export
EAimp <- function(data, weights, outind, reach = "max",
                  transmission.function = "root", power = ncol(data),
                  distance.type = "euclidean", duration = 5, maxl = 5,
                  kdon = 1, monitor = FALSE, threshold = FALSE,
                  deterministic = TRUE, fixedprop = 0) {

  # start computation time
  calc.time <- proc.time()[1]

  # ------- preparation -------

  # number of rows
	n <- nrow(data)

	# number of columns
	p <- ncol(data)

	# set weights to 1 if missing
	if (missing(weights)) {weights <- rep(1, n)}

  # warning if no outliers are passed to function
  if (missing(outind)) {cat("No outlier indicator is given\n")}

	# if there are missing values in outind, set them to FALSE
  if (sum(is.na(outind)) > 0) {
    cat("Missing values in outlier indicator set to FALSE.\n")
    outind[is.na(outind)] <- FALSE
  }

	# find complete and usable (valid information for at least half of var.) records
	response <- !is.na(data)
	complete.records <- apply(response, 1, prod)
	usable.records <- apply(response, 1, sum) >= p/2

	# output to console
	cat("\n Dimensions (n,p):", n, p)
	cat("\n Number of complete records ", sum(complete.records))
	cat("\n Number of records with maximum p/2 variables missing ", sum(usable.records))

	# standardization of weights to sum to sample size
	np <- sum(weights)
	weights <- as.single((n * weights) / np)

	# Calculation of distances (independently of EAdet)
  EA.dist.res <- .EA.dist(data, n = n, p = p, weights = weights, reach = reach,
                          transmission.function = transmission.function,
                          power = power, distance.type = distance.type,
                          maxl = maxl)

  # print progress to console
  if (monitor) {cat("\n\n Distances finished")}

  # check that distances are on same number of obs.
  if (length(EA.dist.res$counterprobs) != n * (n - 1) / 2) {
    cat("\n Distances not on same number of observations")
  }

  # ------- prepare data for imputations -------

  # create new data frame
  imp.data <- data

  # find observations to be imputed
  to.impute <- apply(response, 1, prod) == 0 | outind

  # get indices of observations to be imputed
  imp.ind <- which(to.impute)

  # number of observations to be imputed
  n.imp <- length(imp.ind)
  cat("\n Number of imputands is ", n.imp)

  # reach d0 is set to maximum in order to reach all
  # outliers (max.min.di of .EA.dist)
  d0 <- reach
  cat("\n Reach for imputation is ", d0)

  # ------- start of imputations -------

  # start loop
  for (i in 1:n.imp) {

    # Start epidemic at the observation that needs imputation
    imputand <- imp.ind[i]

    # set imp.time to 1
    imp.time <- 1

    # initialize an infected vector with FALSE for all obs. n
    imp.infected <- rep(FALSE, n)

    # set imputant to infected (TRUE)
    imp.infected[imputand] <- TRUE

    # get indices of missing variables for outlier
    if (outind[imputand]) {
      missvar <- (1:p)
    } else {missvar <- which(!response[imputand, ])}

    # number of missing variables
    n.missvar <- length(missvar)

    # donors must have values for the missing ones and must not be an outlier
    potential.donors <-
      apply(as.matrix(response[ ,missvar]), 1, sum) == n.missvar & !outind

    # imputand itself cannot be a donor (thus, set to FALSE)
    potential.donors[imputand] <- FALSE

    # print progress to console
    if (monitor) {
      cat("\n Imputand: ", imputand)
		  cat(" missvar: ", missvar)
    	cat(" #pot. donors: ", sum(potential.donors))
    }

    # if no donor can impute all missing values of the imputand
    # then maximise number of imputed values
    if (sum(potential.donors) == 0) {
      max.miss.val <- max(apply(as.matrix(response[-imputand, missvar]), 1, sum))
      cat("\n No common donor for all missing values. ", n.missvar - max.miss.val,
          " missing values will remain")
      potential.donors <- apply(response[ ,missvar], 1, sum) == max.miss.val
    }

    # initialize a vector of donors
    donors <- rep(FALSE, n)

    # prepare for loop
    new.imp.infected <- imp.infected
    n.imp.infected <- sum(imp.infected)
    hprod <- rep(1, n)
    imp.infection.time <- rep(0, n)
    imp.infection.time[imp.infected] <- imp.time

    # run epidemic from imputand until at least one potential donor is infected
	  repeat {

	    # exit loop under certain conditions
	    if (sum(donors) >= kdon | imp.time > duration) {break}

	    # add to imp.time
	    imp.time <- imp.time + 1

	    # copy of imp.infected
	    old.imp.infected <- imp.infected

      # what happens here?
	    if(sum(new.imp.infected) > 1) {

	      hprod[!imp.infected] <-
	        hprod[!imp.infected] * apply(sweep(
	          sweep(matrix(EA.dist.res$counterprobs[apply(
	            as.matrix(which(!imp.infected)), 1, .ind.dijs,
	            js = which(new.imp.infected), n = n)], sum(new.imp.infected),
	            n - n.imp.infected), 1, weights[new.imp.infected], "^"), 2,
	          weights[!imp.infected], "^"), 2, prod)

	      } else {

	        if(sum(new.imp.infected) == 1) {

	          hprod[!imp.infected] <-
	            hprod[!imp.infected] * EA.dist.res$counterprobs[apply(
	              as.matrix(which(!imp.infected)), 1, .ind.dijs,
	              js = which(new.imp.infected), n = n)]^(weights[new.imp.infected] *
	                                                       weights[!imp.infected])

	        }
	      }

	    # if deterministic = TRUE, the number of infections is the expected
	    # number and the infected observations are the ones with largest infection
	    # probabilities.
	    if (deterministic) {

	      # at least 1 infection at each step
		    n.to.infect <- max(1, round(sum(1 - hprod[!imp.infected])))

		    # rank is maximum to allow infection
			  imp.infected[!imp.infected] <-
			    rank(1 - hprod[!imp.infected], ties.method = "max") >=
			    n - n.imp.infected - n.to.infect

			  } else {

			    # if TRUE, infect all remaining points with infection probability
			    # above the threshold
			    if (threshold) {
			      imp.infected[!imp.infected] <- hprod[!imp.infected] <= 0.5 ^ (1 / maxl)
			      } else {

			        # if TRUE, a fixed proportion of observations is infected
			        # at each step
			        if (fixedprop > 0) {
			          n.to.infect <- max(1, floor(fixedprop * sum(!imp.infected)))
			          imp.infected[!imp.infected] <-
			            rank(1 - hprod[!imp.infected]) >= n - n.imp.infected - n.to.infect
			        } else {
			          imp.infected[!imp.infected] <-
			            as.logical(rbinom(n - n.imp.infected, 1, 1 - hprod[!imp.infected]))
			        }
			      }
			  }

      # loop results
	    new.imp.infected <- imp.infected & (!old.imp.infected)
		  n.imp.infected <- sum(imp.infected)
		  imp.infection.time[new.imp.infected] <- imp.time
      donors <- potential.donors & imp.infected

      # print progress to console
      if (monitor) {cat(", donors: ", which(donors))}

      # start next evaluation of loop
		  next
	  }

    # if no donors are found, print message to console
    if (imp.time > duration) {

      cat("No donor found")

      } else {

        # if there are more than one donor, sample from donors
        if (sum(donors) > 1) {

          don.imp <- sample(which(donors), 1)

          } else {

            don.imp <- which(donors)

          }

        # print progress to console
        if (monitor) {cat(". #donors: ", sum(donors), ", chosen donor: ", don.imp)}

        # IMPUTE data with data from donor
        imp.data[imputand, missvar] <- data[don.imp, missvar]

      }
  }

  # stop computation time
  calc.time <- round(proc.time()[1] - calc.time, 5)

  # output to console the remaining number of missing values
  cat("\n\n Number of remaining missing values is ", sum(is.na(imp.data)))

  # ------- results -------

  # prepare output
	EAimp.r <- list(sample.size = n, number.of.variables = p,
	                n.complete.records = sum(complete.records),
	                n.usable.records = sum(usable.records),
	                duration = duration, reach = d0, threshold = threshold,
	                deterministic = deterministic, computation.time = calc.time)

  # return output
  return(invisible(list(output = EAimp.r, imputed.data = imp.data)))

  # output to console
	cat("\n", "EA imputation has finished in", calc.time, "seconds.", "\n")
}

