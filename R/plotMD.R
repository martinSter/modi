#' QQ-Plot of Mahalanobis distances
#'
#' QQ-plot of (squared) Mahalanobis distances vs. scaled F-distribution (or a scaled chisquare distribution).
#' In addition, two default cutpoints are proposed.
#'
#' Scaling of the F-distribution as \code{median(dist)*qf((1:n)/(n+1), p, n-p)/qf(0.5, p, n-p)}.
#' First default cutpoint is \code{median(dist)*qf(alpha, p, n-p)/qf(0.5, p, n-p)} and the second default
#' cutpoint is the alpha quantile of the Mahalanobis distances.
#'
#' @param dist a vector of Mahalanobis distances.
#' @param p the number of variables involved in the Mahalanobis distances.
#' @param alpha a probability for cut-off, usually close to 1.
#' @param chisquare a logical indicating the the chisquare distribution should be used instead of the F-distribution.
#' @return \item{hmed}{first proposed cutpoint based on F-distribution}
#' @return \item{halpha}{second proposed cutpoint (alpha-quantile)}
#' @return \item{QQ-plot}{}
#' @author Beat Hulliger
#' @references Little, R. & Smith, P. (1987) Editing and imputation for quantitative survey data,
#' Journal of the American Statistical Association, 82, 58-68
#' @examples
#' data(bushfirem, bushfire.weights)
#' det.res <- TRC(bushfirem, weights = bushfire.weights)
#' PlotMD(det.res$dist, ncol(bushfirem))
#' @export
#' @importFrom stats qf qchisq median
#' @importFrom graphics plot abline title
PlotMD <- function(dist, p, alpha = 0.95, chisquare = FALSE) {

  # QQ-Plot of (squared) Mahalanobis Distance
	# dist[is.na(dist)]<-max(dist,na.rm=TRUE) # set missing distances to max

  # remove missing values from dist
  dist <- dist[!is.na(dist)]

  # save length of dist
	n <- length(dist)

	# compute quantiles
	if (!chisquare) {
	  # F-distribution
    x <- median(dist) * qf((1:n) / (n + 1), p, n - p) / qf(0.5, p , n - p)
	} else {
	  # Chisquare distribution
	  x <- median(dist) * qchisq((1:n) / (n + 1), p) / qchisq(0.5, p)
	}

	# sort dist in increasing order
  y <- sort(dist)

  # create plot and comput cutpoints
  if (!chisquare) {
    # F-distribution
    plot(x, y, ylab = "MD-quantiles", xlab = "F-quantiles")
    hmed <- median(dist) * qf(alpha, p, n - p) / qf(0.5, p, n - p)
    halpha <- sort(dist)[floor(alpha * n)]
  } else {
    # Chisquare distribution
    plot(x, y, ylab = "MD-quantiles", xlab = "Chisquare-quantiles")
    hmed <- median(dist) * qchisq(alpha, p) / qchisq(0.5, p)
    halpha <- sort(dist)[floor(alpha * n)]
  }

  # abline hmed and halpha
  abline(h = hmed, lty = 1)
  abline(h = halpha, lty = 2)

  # set titles
  title(main = "QQ-Plot of Mahalanobis distances")
  title(sub = paste("alpha = ", alpha, ", hmed = ", round(hmed, 3),
                    ", halpha = ", round(halpha, 3),
                    " n.miss.dist = ", sum(is.na(dist))))

  # return two cutpoints
  return(list(hmed = hmed, halpha = halpha))
}
