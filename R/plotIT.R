#' Plot of  infection times of the EA algorithm
#'
#' The (weighted) cdf of infection times is plotted. The infection times
#' vertical jumps of the cdf is shown by the points with the
#' same infection times stacked vertically and respecting the weights.
#'
#' The infection times of \code{EAdet} are the main input. In addition the weights
#' may be needed. The default cutpoint from \code{EAdet} may be used for the cutpoint.
#' Points that are never infected have a missing infection time. These missing infection times
#' are (temporarily) imputed by 1.2 times the maximum infection time
#' to show them on the plot marked with an x.
#'
#' @param infection.time vector of infection.times of the observations
#' @param weights vector of (survey) weights of the observations
#' @param cutpoint a cutpoint to for declaring outliers
#'
#' @author Beat Hulliger
#' @examples
#' it <- c(rep(NA, 3), rep(1:7, times=c(1, 4, 10, 8, 5, 3, 2)))
#' wt <- rep(c(1,2,5), times=12)
#' plotIT(it, wt, 6)
#' @export
#' @importFrom graphics plot abline title points
plotIT <- function(infection.time, weights, cutpoint)
{
  weights <- weights/sum(weights)
  plot.time <-infection.time
  sel <- is.na(infection.time)
  plot.time[is.na(infection.time)] <- ceiling(1.2 * max(infection.time, na.rm=TRUE))
  ord <- order(plot.time)
  x <- plot.time[ord]
  y <- cumsum(weights[ord])
  plot(y~x, xlab="infection time",  ylab = "(weighted) cdf of infection time")
  points(x[sel[ord]],y[sel[ord]], col="red")
  abline(v=cutpoint)
}
