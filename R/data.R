#' Bushfire scars.
#'
#' The bushfire data set was used by Campbell (1984, 1989) to locate bushfire scars.
#' The dataset contains satelite measurements on five frequency bands, corresponding
#' to each of 38 pixels.
#'
#' The data contains an outlying cluster of observations 33 to 38 a second outlier
#' cluster of observations 7 to 11 and a few more isolated outliers, namely observations
#' 12, 13, 31 and 32.
#'
#' For testing purposes weights are provided:
#' \code{bushfire.weights <- rep(c(1,2,5), length = nrow(bushfire))}
#'
#' @format A data frame with 38 rows and 5 variables.
#' @references Campbell, N. (1989) Bushfire Mapping using NOAA AVHRR Data. Technical Report.
#' Commonwealth Scientific and Industrial Research Organisation, North Ryde.
#' @examples
#' data(bushfire)
"bushfire"

#' Bushfire scars with missing data.
#'
#' The bushfire data set was used by Campbell (1984, 1989) to locate bushfire scars.
#' The dataset contains satelite measurements on five frequency bands, corresponding
#' to each of 38 pixels. However, this dataset contains missing values.
#'
#' The data contains an outlying cluster of observations 33 to 38 a second outlier
#' cluster of observations 7 to 11 and a few more isolated outliers, namely observations
#' 12, 13, 31 and 32.
#'
#' \code{bushfirem} is created from bushfire by setting a proportion of 0.2 of the values
#' to missing.
#'
#' For testing purposes weights are provided:
#' \code{bushfire.weights <- rep(c(1,2,5), length = nrow(bushfire))}
#'
#' @format A data frame with 38 rows and 5 variables.
#' @references Campbell, N. (1989) Bushfire Mapping using NOAA AVHRR Data. Technical Report.
#' Commonwealth Scientific and Industrial Research Organisation, North Ryde.
#' @examples
#' data(bushfirem)
"bushfirem"

#' Weights for Bushfire scars.
#'
#' The bushfire data set was used by Campbell (1984, 1989) to locate bushfire scars.
#' The dataset contains satelite measurements on five frequency bands, corresponding
#' to each of 38 pixels.
#'
#' For testing purposes, \code{bushfire.weights} provides artificial weights created
#' according to: \code{bushfire.weights <- rep(c(1,2,5), length = nrow(bushfire))}
#'
#' @format A vector of length 38.
#' @references Campbell, N. (1989) Bushfire Mapping using NOAA AVHRR Data. Technical Report.
#' Commonwealth Scientific and Industrial Research Organisation, North Ryde.
#' @examples
#' data(bushfire.weights)
"bushfire.weights"

#' Sample Environment Protection Expenditure Survey.
#'
#' The sepe data set is a sample of the pilot survey in 1993 of the Swiss Federal Statistical
#' Office on environment protection expenditures of Swiss private economy in the previous
#' accounting year. The units are enterprises, the monetary variables are in thousand Swiss
#' Francs (CHF). From the original sample a random subsample was chosen of which certain
#' enterprises were excluded for confidentiality reasons. In addition, noise has been added
#' to certain variables, and certain categories have been collapsed. The data set has missing
#' values. The data set has first been prepared for the EU FP5 project EUREDIT and later been
#' data protected for educational purposes.
#'
#' The sample design is stratified random sampling with different sampling rates. Use package
#' survey or sampling to obtain correct point and variance estimates. In addition a ratio
#' estimator may be built using the variable popemple which gives the total employment per
#' activity.
#'
#' There are two balance rules: the subtotals of the investment variables should
#' sum to totinvto and the expenditure subtotals should sum to totexpto.
#'
#' The missing values stem from the survey itself. In the actual survey the missing
#' values were declared as 'guessed' rather than copied from records.
#'
#' The sampling weight weight is adjusted for non-response in the stratum,
#' i.e. \code{weight=popsize/sampsize}.
#'
#' @format A data frame with 675 rows and 23 variables:
#' \describe{
#'   \item{idnr}{identifier (anonymous)}
#'   \item{exp}{categorical variable where 1 = 'non-zero total expenditure' and
#'   2 = 'zero total expenditure, and 3 = 'no answer'}
#'   \item{totinvwp}{total investment for water protection}
#'   \item{totinvwm}{total investment for waste management}
#'   \item{totinvap}{total investment for air protection}
#'   \item{totinvnp}{total investment for noise protection}
#'   \item{totinvot}{total investement for other environmental protection}
#'   \item{totinvto}{overall total investment in all environmental protection areas}
#'   \item{totexpwp}{total current expenditure in environmental protectiona area water protection}
#'   \item{totexpwm}{total current expenditure in environmental protectiona area waste management}
#'   \item{totexpap}{total current expenditure in environmental protectiona area air protection}
#'   \item{totexpnp}{total current expenditure in environmental protectiona area noise protection}
#'   \item{totexpot}{total current expenditure in other environmental protection}
#'   \item{totexpto}{overall total current expenditure in all environmental protection}
#'   \item{subtot}{total subsidies for environmental protection received}
#'   \item{rectot}{total receipts from environmental protection}
#'   \item{employ}{number of employees}
#'   \item{sizeclass}{size class (according to number of employees)}
#'   \item{stratum}{stratum number of sample design}
#'   \item{activity}{code of economic activity (aggregated)}
#'   \item{popsize}{number of enterprises in the population-stratum}
#'   \item{popempl}{number of employees in population activity group}
#'   \item{weight}{sampling weight (for extrapolation to the population)}
#' }
#' @references
#' Swiss Federal Statistical Office (1996), Umweltausgaben und -investitionen in der
#' Schweiz 1992/1993, Ergebnisse einer Pilotstudie.
#'
#' Charlton, J. (ed.), Towards Effective Statistical Editing and Imputation Strategies -
#' Findings of the Euredit project, unpublished manuscript available from Eurostat
#' and \url{http://www.cs.york.ac.uk/euredit/}.
#' @examples
#' data(sepe)
"sepe"
