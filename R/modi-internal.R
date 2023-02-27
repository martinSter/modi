#' Utility function for EAdet and EAimp
#'
#' Calculation of distances for EPIDEMIC Algorithm for
#' multivariate outlier detection and imputation
#'
#' @param data a data frame or matrix with data.
#' @param n number of rows.
#' @param p number of columns.
#' @param weights a vector of positive sampling weights.
#' @param reach if \code{reach = "max"} the maximal nearest
#' neighbor distance is used as the basis for the transmission
#' function, otherwise the weighted
#' @param transmission.function form of the transmission function
#' of distance d: \code{"step"} is a heaviside function which jumps
#' to \code{1} at \code{d0}, \code{"linear"} is linear between
#' \code{0} and \code{d0}, \code{"power"} is \code{(beta*d+1)^(-p)}
#' for \code{p = ncol(data)} as default, \code{"root"} is the function
#' \code{1-(1-d/d0)^(1/maxl)}.
#' @param power sets \code{p = power}.
#' @param distance.type distance type in function \code{dist()}.
#' @param maxl maximum number of steps without infection.
#' @keywords internal
#' @author Beat Hulliger
#' @export
#' @importFrom stats dist
EA.dist <- function(data, n, p, weights, reach,
                    transmission.function, power,
                    distance.type, maxl) {

    # save distances, i.e. the later counterprobabilities
    counterprobs <- as.single(dist(data, method = distance.type))

	# The dist function handles missing values correctly except
	# if there is no overlap (see counterprob)
	min.di <- rep(0, n)
	means.di <- rep(0, n)

	# Will be used for the sample spatial median and for d0
	for (i in 1:n) {
	    di <- counterprobs[ind.dijs(i, 1:n, n)]
		  min.di[i] <- nz.min(di)

		  # weighted mean of distances to account for missing distances
		  means.di[i] <- sum(di * weights[-i], na.rm = TRUE) /
            sum(weights[-i][!is.na(di)])
	}

	# Sample spatial median
    # Restrict to usable observations (for sample spatial median)
	usable.records <- apply(!is.na(data), 1, sum) >= p / 2
	means.di.complete <- means.di
	means.di.complete[!usable.records] <- NA
	sample.spatial.median.index <- which(means.di.complete ==
        min(means.di.complete, na.rm = TRUE))[1]

	# Determine tuning distance d0
	max.min.di <- max(min.di, na.rm = TRUE)
    d0 <- if (is.numeric(reach)) {
        reach
    } else {
        switch(reach,
            "max" = min(max.min.di, 2 * sqrt(p)),
            "quant" = min(weighted.quantile(min.di, w = weights,
                prob = 1 - (p + 1) / n), 2 * sqrt(p)),
            stop("argument 'reach' is not defined\n"))
    }

	# Calculation of counterprobabilities
    # counterprobabilities stocked in distances vector to save memory
    # counterprobabilities are set to 1 if missing
	if (n%%2 == 0) {
	    l.batch <- n - 1
	    n.loops <- n / 2
	} else {
	    l.batch <- n
		n.loops <- (n - 1) / 2
	}

	if (transmission.function == "step") {
        for (i in 1:n.loops) {
	        dij <- counterprobs[(i - 1) * l.batch + (1:l.batch)]
	        dij <- as.single(dij > d0)
	        dij[is.na(dij)] <- 1
	        counterprobs[(i - 1) * l.batch + (1:l.batch)] <- as.single(dij)
	    }

    } else {
        if(transmission.function == "linear") {
			for(i in 1:n.loops) {
			  dij <- counterprobs[(i - 1) * l.batch + (1:l.batch)]
			  dij <- 1 - pmax(0, d0 - dij) / d0
			  dij[is.na(dij)] <- 1
			  counterprobs[(i - 1) * l.batch + (1:l.batch)] <- as.single(dij)
			}

	    } else {
            if (transmission.function == "power") {
                beta <- as.single((0.01^(-1 / power) - 1) / d0)

                for (i in 1:n.loops) {
                    dij <- counterprobs[(i - 1) * l.batch + (1:l.batch)]
                    dij <- 1 - (beta * dij + 1)^(-power)
                    dij <- ifelse(dij > d0, 1, dij)
                    dij[is.na(dij)] <- 1
                    counterprobs[(i - 1) * l.batch + (1:l.batch)] <- as.single(dij)
                  }

            } else { # default transmission function is the root function
                for (i in 1:n.loops) {
                    dij <- counterprobs[(i - 1) * l.batch + (1:l.batch)]
                    dij <- 1 - (1 - dij / d0)^(1 / power)
                    dij <- ifelse(dij > d0, 1, dij)
                    dij[is.na(dij)] <- 1
                    counterprobs[(i - 1) * l.batch + (1:l.batch)] <- as.single(dij)
                }
            }
        }
    }

	return(invisible(list(
	  sample.spatial.median.index = sample.spatial.median.index,
	  max.min.di = max.min.di,
	  transmission.distance = d0,
	  min.di = min.di,
	  counterprobs = counterprobs)))
}



#' Non-zero non-missing minimum function
#'
#' Returns the non-zero non-missing minimum.
#'
#' @param x vector of data.
#' @keywords internal
#' @author Beat Hulliger
#' @export
nz.min <- function(x) {
  if (sum(!is.na(x)) == 0) {
    nz.min.temp <- NA
  } else {
    nz.min.temp <- min(x[x != 0], na.rm = TRUE)
  }
  return(nz.min.temp)
}



#' Addressing function for Epidemic Algorithm
#'
#' Utility function for Epidemic Algorithm.
#'
#' @param i index i.
#' @param j index j.
#' @param n number of rows.
#' @keywords internal
#' @author Cédric Béguin
#' @export
ind.dij <- function(i, j, n) {
  return((i - 1) * n - ((i + 1) * i) / 2 + j)
}



#' Addressing function for Epidemic Algorithm
#'
#' Utility function for Epidemic Algorithm.
#'
#' @param i index i.
#' @param js indexes js.
#' @param n number of rows.
#' @keywords internal
#' @author Cédric Béguin
#' @export
ind.dijs <- function(i, js, n){
  indices <- c(ind.dij(js[js < i], i, n), ind.dij(i, js[js > i], n))
  return(indices[!is.na(indices)])
}



#' Utility function for TRC.R among others
#'
#' Sum of weights for observations < value (if lt = TRUE)
#' or observations=value (if lt = FALSE).
#'
#' @param observations vector of observations.
#' @param weights vector of weights.
#' @param value value.
#' @param lt either \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @author Beat Hulliger
#' @export
weightsum <- function(observations, weights, value, lt = TRUE) {
  if (lt) {
    return(sum(weights * (observations < value), na.rm = TRUE))
  } else {
    return(sum(weights * (observations == value), na.rm = TRUE))
  }
}



#' Sweep operator
#'
#' Definition of the sweep and reverse-sweep operator (Schafer pp 159-160)
#'
#' @param M a matrix.
#' @param k column.
#' @param reverse either \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @author Beat Hulliger
#' @export
sweep.operator <- function(M, k, reverse = FALSE) {
  if (reverse) {
    Gjk <- M[ ,k]
    Hkk <- -1 / M[k,k]
    M <- M + (Gjk %*% t(Gjk)) * Hkk
    M[k, ] <- M[ ,k] <- Gjk * Hkk
    M[k,k] <- Hkk
    return(M)
  } else {
    Gjk <- M[ ,k]
    Hkk <- 1 / M[k,k]
    M <- M - (Gjk %*% t(Gjk)) * Hkk
    M[k, ] <- M[ ,k] <- Gjk * Hkk
    M[k,k] <- -Hkk
    return(M)
  }
}



#' psi-function
#'
#' Defined in Little and Smith for ER algorithm
#'
#' @param d vector of distances.
#' @param present present.
#' @param psi.par parameters.
#' @keywords internal
#' @author Beat Hulliger
#' @export
psi.lismi <- function(d, present, psi.par = c(2, 1.25)) {
  a <- psi.par[1]
  b <- psi.par[2]
  d0 <- sqrt(present) + a / 2
  return(ifelse(d - d0 <= 0, d, d0 * exp(-(d - d0)^2 / (2 * b^2))))
}



#' EM for multivariate normal data
#'
#' This version of EM does not contain the computation of the
#' observed sufficient statistics, they will be computed in the
#' main program of BEM and passed as parameters as well as the
#' statistics on the missingness patterns.
#'
#' @param data matrix or dataframe with data.
#' @param weights vector of weights.
#' @param n number of rows.
#' @param p number of columns.
#' @param s.counts s.counts.
#' @param s.id s.id.
#' @param S S.
#' @param T.obs T.obs.
#' @param start.mean initial center.
#' @param start.var initial variance.
#' @param numb.it numb.it.
#' @param Estep.output Estep.output.
#' @keywords internal
#' @author Beat Hulliger
#' @export
EM.normal <- function(data, weights = rep(1, nrow(data)), n = sum(weights),
                       p = ncol(data), s.counts, s.id, S, T.obs,
                       start.mean = rep(0, p), start.var = diag(1, p),
                       numb.it = 10, Estep.output = FALSE) {

  # ------- Initialization -------

  # Creates theta which is the matrix form of the initial
  # parameter used by EM
  theta <- matrix(0, p + 1, p + 1)
  theta[1,1] <- -1
  theta[1,2:(p + 1)] <- theta[2:(p + 1), 1] <- start.mean
  theta[2:(p + 1),2:(p + 1)] <- start.var

  # ------- Iterations of EM -------
  for (boucle in 1:numb.it) {

    if (Estep.output) {cat("E-step ", boucle, "\n")}

    # ------- The E-step -------

    # Initializing T.tot to T.obs
    T.tot <- T.obs

    # Start loop on missing patterns s from 1 to S
    for (s in 1:S) {

      # Identification of the indices of x.mis and x.obs
      x.mis.id <- (1:p)[is.na(data[s.id[s], ])]
      x.obs.id <- (1:p)[-x.mis.id]

      # Sweep of theta over the indices of x.obs
      C.s <- theta

      for (k in x.obs.id) {
        if (C.s[k + 1, k + 1] != 0) {
          C.s <- sweep.operator(C.s, k + 1)
        }
      }

      # Start loop over all observations x having missing pattern s
      for (i in 1:s.counts[s]) {
        if (s == 1) {
          x <- data[i, ]
          weight <- weights[i]
        } else {
          x <- data[s.id[s - 1] + i, ]
          weight <- weights[s.id[s - 1] + i]
        }

        # Computation of x.star = E(x.mis|x.obs)
        x.star <- x

        for (k in x.mis.id) {
          x.star[k] <- C.s[1, k + 1] + sum(C.s[x.obs.id + 1, k + 1] * x[x.obs.id])
        }

        # Updating T.tot
        T.tot[1, ] <- T.tot[ ,1] <- T.tot[1, ] + weight * c(1, x.star)
        T.tot[2:(p + 1), 2:(p + 1)] <- T.tot[2:(p + 1), 2:(p + 1)] + weight * (x.star %*% t(x.star))
        T.tot[x.mis.id + 1, x.mis.id + 1] <- T.tot[x.mis.id + 1, x.mis.id + 1] + weight * C.s[x.mis.id + 1, x.mis.id + 1]

      }

    }

    # ------- The M-step -------
    theta <- sweep.operator(T.tot / n, 1)

  }

  # ------- End of EM -------
  return(theta)

}



#' Utility for ER function
#'
#' The \code{ER} function is an implementation of the
#' ER-algorithm of Little and Smith (1987).
#'
#' @param data matrix or dataframe with data.
#' @param weights vector of weights.
#' @param psi.par parameters for psi function.
#' @param np np.
#' @param p number of columns.
#' @param s.counts s.counts.
#' @param s.id s.id.
#' @param S S.
#' @param missing.items missing items.
#' @param nb.missing.items number of missing items.
#' @param start.mean initial center.
#' @param start.var initial variance.
#' @param numb.it number of iterations.
#' @param Estep.output Estep.output.
#' @param tolerance tolerance.
#' @keywords internal
#' @author Beat Hulliger
#' @export
ER.normal <- function(data, weights = rep(1, nrow(data)),
                       psi.par = c(2, 1.25), np = sum(weights),
                       p = ncol(data), s.counts, s.id, S, missing.items,
                       nb.missing.items, start.mean = rep(0, p),
                       start.var = diag(1, p), numb.it = 10,
                       Estep.output = FALSE, tolerance = 1e-06) {

  # ------- Initialization -------

  # Creates theta which is the matrix form of the
  # initial parameter used by EM
  theta <- matrix(0, p + 1, p + 1)
  theta[1,1] <- -1
  theta[1, 2:(p + 1)] <- theta[2:(p + 1), 1] <- start.mean
  theta[2:(p + 1), 2:(p + 1)] <- start.var

  if (Estep.output) {cat("\n", "theta: \n", theta)}

  break.flag <- FALSE

  # Initialisation of robustness weights to 1
  rob.weights <- rep(1, nrow(data))
  np.hat <- np

  # ------- Iterations of EM -------

  for (boucle in 1:numb.it) {

    if (Estep.output) {cat("\n E-step ", boucle)}

    # ------- The E-step -------

    # Initializing T.tot and storing of old.theta
    T.tot <- matrix(0, p + 1, p + 1)
    old.theta <- theta

    # Start loop on missing patterns s from 1 to S
    for (s in 1:S) {

      # Identification of the indices of x.mis and x.obs
      x.mis.id <- (1:p)[is.na(data[s.id[s], ])]
      x.obs.id <- (1:p)[-x.mis.id]

      # Sweep of theta over the indices of x.obs
      C.s <- theta

      for (k in x.obs.id) {
        if (C.s[k + 1, k + 1] != 0) {
          C.s <- sweep.operator(C.s, k + 1)
        }
      }

      # Start loop over all observations x having missing pattern s
      for (i in 1:s.counts[s]) {
        if (s == 1) {
          x <- data[i, ]
          weight <- weights[i]
          rob.weight <- rob.weights[i]
        } else {
          x <- data[s.id[s - 1] + i, ]
          weight <- weights[s.id[s - 1] + i]
          rob.weight <- rob.weights[s.id[s - 1] + i]
        }

        # Computation of x.star = E(x.mis|x.obs)
        x.star <- x

        for (k in x.mis.id) {
          x.star[k] <- C.s[1, k + 1] + sum(C.s[x.obs.id + 1, k + 1] * x[x.obs.id])
        }

        # robustification with robustness weights from last iteration
        x.star <- x.star

        # Updating T.tot
        T.tot[1, ] <- T.tot[ ,1] <- T.tot[1, ] + rob.weight * weight * c(1, x.star)
        T.tot[2:(p + 1), 2:(p + 1)] <- T.tot[2:(p + 1), 2:(p + 1)] + rob.weight * weight * (x.star %*% t(x.star))
        T.tot[x.mis.id + 1, x.mis.id + 1] <- T.tot[x.mis.id + 1, x.mis.id + 1] + rob.weight * weight * C.s[x.mis.id + 1, x.mis.id + 1]
      }
    }

    # ------- The M-step -------

    np.hat <- sum(weights * rob.weights)
    theta <- sweep.operator(T.tot / np.hat, 1)

    if (Estep.output) {cat("\n", "theta: \n", theta, "\n")}

    # Computation of Mahalanobis distances (marginal version)
    ER.mean <- theta[1, 2:(p + 1)]
    ER.var <- theta[2:(p + 1), 2:(p + 1)]
    indices <- (!missing.items[1, ])
    dist <- mahalanobis(data[1:s.id[1], indices, drop = FALSE],
                        ER.mean[indices],
                        ER.var[indices, indices]) * p / (p - nb.missing.items[1])

    if (S > 1) {

      for (i in 2:S) {
        indices <- (!missing.items[i, ])
        dist <-
          c(dist,
            mahalanobis(data[(s.id[i - 1] + 1):s.id[i], indices, drop = FALSE],
                        ER.mean[indices],
                        ER.var[indices, indices, drop = FALSE]) * p / (p - nb.missing.items[i]))
      }

    }

    # new robustness weights
    rob.weights <-
      ifelse(dist > 0, psi.lismi(sqrt(dist), apply(!is.na(data), 1, sum), psi.par = psi.par) / sqrt(dist), 1)

    if (Estep.output) {cat("\n Delta: ", signif(max(abs(theta-old.theta)),8),"\n")}

    if (max(abs(theta - old.theta)) < tolerance) {
      break.flag <- T; break
    }

  }

  return(list(theta = theta, dist = dist,
              rob.weights = rob.weights, convergence = break.flag))

}
