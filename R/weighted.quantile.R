#' Quantiles of a weighted cdf
#'
#' A weighted cdf is calculated and quantiles are evaluated. Missing values are discarded.
#'
#' Weighted linear interpolation in case of non-unique inverse. Gives a warning when the
#' contribution of the weight of the smallest observation to the total weight is larger
#' than \code{prob}.
#'
#' @param x a vector of data.
#' @param w a vector of (sampling) weights.
#' @param prob the probability for the quantile.
#' @param plot if \code{TRUE}, the weighted cdf is plotted.
#' @return The quantile according to \code{prob} (by default it returns the weighted median).
#' @author Beat Hulliger
#' @note No variance calculation.
#' @seealso \href{https://www.rdocumentation.org/packages/survey/versions/3.33-2/topics/svyquantile}{svyquantile}
#' @examples
#' x <- rnorm(100)
#' x[sample(1:100, 20)] <- NA
#' w <- rchisq(100, 2)
#' weighted.quantile(x, w, 0.2, TRUE)
#' @export
#' @importFrom stats quantile
#' @importFrom graphics plot
weighted.quantile <- function(x, w, prob = 0.5, plot = FALSE) {

  # if weights are missing, we return simple quantile function
  # otherwise, we remove the weights where the data are missing
  if(missing(w)) {

    return(quantile(x, prob, na.rm = TRUE))

  } else {

    w <- w[!is.na(x)]

  }

  # remove missing data
  x <- x[!is.na(x)]

  # sort the data and the weights in ascending order
	ord <- order(x)
	w <- w[ord]
	x <- x[ord]

  # compute weighted cdf
	w.ord <- cumsum(w) / sum(w)

	# get index with length of data vector x
	index <- seq_along(x)

	# plot weighted cdf (w.ord) if plot = TRUE
	if (plot) {plot(x, w.ord, type = "s")}

	# if the min of the cdf is larger than prob,
	# then the quantile is based on only one observation (Warning)
	# thus, set index of lower quant to 1 (first observation)
	# otherwise, take max index where w.ord smaller or equal to prob
	if(min(w.ord) > prob) {

	  warning('Dominance of one observation!')
		lower.k.quant <- 1

	} else {

	  lower.k.quant <- max(index[w.ord <= prob])

	}

	# set index of upper quant to min index where w.ord greater than prob
	upper.k.quant <- min(index[w.ord > prob])

	# if inverse is non-unique, return weighted linear interpolation
	if(w.ord[lower.k.quant] < prob) {

	  # returns quantile according to upper.k.quant
	  return(x[upper.k.quant])

	} else {

	  # returns weighted linear interpolation
	  return((w[lower.k.quant] * x[lower.k.quant] + w[upper.k.quant] * x[upper.k.quant]) /
	           (w[lower.k.quant] + w[upper.k.quant]))

	}
}

