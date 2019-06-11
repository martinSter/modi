#' Winsorization followed by imputation
#'
#' Winsorization of outliers according to the Mahalanobis distance
#' followed by an imputation under the multivariate normal model.
#' Only the outliers are winsorized. The Mahalanobis distance MDmiss
#' allows for missing values.
#'
#' It is assumed that \code{center}, \code{scatter} and \code{outind}
#' stem from a multivariate outlier detection algorithm which produces
#' robust estimates and which declares outliers observations with a large
#' Mahalanobis distance. The cutpoint is calculated as the least (unsquared)
#' Mahalanobis distance among the outliers. The winsorization reduces the
#' weight of the outliers:
#' \deqn{\hat{y}_i = \mu_R + (y_i - \mu_R) \cdot c/d_i}
#' where \eqn{\mu_R} is the robust center and \eqn{d_i} is the (unsquared) Mahalanobis
#' distance of observation i.
#'
#' @param data a data frame with the data.
#' @param center (robust) estimate of the center (location) of the observations.
#' @param scatter (robust) estimate of the scatter (covariance-matrix) of the
#' observations.
#' @param outind logical vector indicating outliers with 1 or TRUE for outliers.
#' @param seed seed for random number generator.
#' @return \code{Winsimp} returns a list whose first component \code{output} is a
#' sublist with the following components:
#' \describe{
#'   \item{\code{cutpoint}}{Cutpoint for outliers}
#'   \item{\code{proc.time}}{Processing time}
#'   \item{\code{n.missing.before}}{Number of missing values before imputation}
#'   \item{\code{n.missing.after}}{Number of missing values after imputation}
#' }
#' The further component returned by \code{winsimp} is:
#' \describe{
#'   \item{\code{imputed.data}}{Imputed data set}
#' }
#' @author Beat Hulliger
#' @references Hulliger, B. (2007), Multivariate Outlier Detection and Treatment
#' in Business Surveys, Proceedings of the III International Conference on
#' Establishment Surveys, Montréal.
#' @seealso \code{\link{MDmiss}}. Uses \code{\link[norm]{imp.norm}}.
#' @examples
#' data(bushfirem, bushfire.weights)
#' det.res <- TRC(bushfirem, weight = bushfire.weights)
#' imp.res <- Winsimp(bushfirem, det.res$output$center, det.res$output$scatter, det.res$outind)
#' print(imp.res$output)
#' @export
Winsimp <- function(data, center, scatter, outind, seed = 1000003) {

  # start computation time
  calc.time <- proc.time()

  # convert outind to logical
  outind <- as.logical(outind)

  # transform data to matrix
  data.wins <- as.matrix(data)

  # compute Mahalanobis distance (not squared)
  MD <- sqrt(MDmiss(data.wins, center, scatter))

  # define minimal MD of outliers as cutpoint
  cutpoint <- min(MD[outind], na.rm = TRUE)

  # robustness weight (the larger distance of MD to cutpoint, the smaller weight)
  u <- ifelse(MD <= cutpoint, 1, cutpoint/MD)

  ###  winsorization for outliers (only outliers!)
  # first sweep subtracts mean from outliers
  # secord sweep multiplies robustness weights with outliers
  # third sweep adds mean again
  data.wins[outind, ] <- as.matrix(
    sweep(sweep(sweep(data[outind, ], 2, center, '-'),
                1, u[outind], '*'), 2, center, '+'))

  ### imputation for missing values with norm package

  # set seed
  norm::rngseed(seed)

  # sorts rows of data.wins by missingness patterns, and centers/scales columns
  s <- norm::prelim.norm(data.wins)

  # impute missing multivariate normal data
  data.imp <- norm::imp.norm(s, norm::makeparam.norm(s,list(center, scatter)), data.wins)

  # print to console if there are missing values
  if (sum(is.na(data.imp)) > 0) {
    cat("There are missing values in the imputed data set.\n")
  }

  # stop computation time
	calc.time <- proc.time() - calc.time

	# prepare output
  res <- list(cutpoint = cutpoint,
              proc.time = calc.time,
              n.missing.before = sum(is.na(data)),
              n.missing.after = sum(is.na(data.imp)),
              imputed.data = data.imp)

  class(res) <- "Winsimp.r"
  res
}

