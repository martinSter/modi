#' Weighted univariate variance coping with missing values
#'
#' This function is analogous to weighted.mean.
#'
#' The weights are standardised such that \eqn{\sum_{observed} w_i} equals the number of observed
#' values in \eqn{x}. The function calculates \deqn{\sum_{observed} w_i(x_i -
#' weighted.mean(x, w, na.rm = TRUE))^2/((\sum_{observed} w_i) - 1)}
#'
#' @param x a vector with data.
#' @param w a vector of positive weights (may not have missings where x is observed).
#' @param na.rm if \code{TRUE} remove missing values.
#' @return The weighted variance of \code{x} with weights \code{w} (with missing values removed
#' when \code{na.rm = TRUE}).
#' @author Beat Hulliger
#' @seealso \href{http://stat.ethz.ch/R-manual/R-devel/library/stats/html/weighted.mean.html}{weighted.mean}
#' @examples
#' x <- rnorm(100)
#' x[sample(1:100, 20)] <- NA
#' w <- rchisq(100, 2)
#' weighted.var(x, w, na.rm = TRUE)
#' @export
weighted.var <- function (x, w, na.rm = FALSE) {

  if (missing(w)) {
    # if weights are missing, we set them to 1
    w <- rep.int(1, length(x))
  } else if (length(w) != length(x)) {
    # if data and weights do not have same length, throw error
    stop("x and w must have the same length")
  }

  # if weights are negative, throw error
  if (min(w)<0) stop("there are negative weights")

  # if weights are integers, define as numeric
  if (is.integer(w)) w <- as.numeric(w)

  # if na.rm = TRUE, remove missing values
  if (na.rm) {
    w <- w[obs.ind <- !is.na(x)]
    x <- x[obs.ind]
  }

  # standardize weights such that the sum of weights equals the number of observed values
  w <- w * length(w) / sum(w)

  return(sum(w * (x - weighted.mean(x, w))^2) / (sum(w) - 1))
}

