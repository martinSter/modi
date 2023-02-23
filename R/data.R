#' Bushfire scars.
#'
#' The bushfire data set was used by Campbell (1984, 1989) to locate bushfire scars.
#' The dataset contains satellite measurements on five frequency bands, corresponding
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
#' The dataset contains satellite measurements on five frequency bands, corresponding
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
#' The dataset contains satellite measurements on five frequency bands, corresponding
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
#'   \item{totinvot}{total investment for other environmental protection}
#'   \item{totinvto}{overall total investment in all environmental protection areas}
#'   \item{totexpwp}{total current expenditure in environmental protection area water protection}
#'   \item{totexpwm}{total current expenditure in environmental protection area waste management}
#'   \item{totexpap}{total current expenditure in environmental protection area air protection}
#'   \item{totexpnp}{total current expenditure in environmental protection area noise protection}
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

#' Living Standards Measurement Survey Albania 2012
#'
#' The dataset is an extended version of the public micro data file of the LSMS 2012 of
#' Albania (\url{http://www.instat.gov.al/media/1547/lsms_2012_eng.rar}, accessed 14 August 2020).
#' Documentation of the LSMS 2012 of Albania is from the
#' World Bank (\url{https://microdata.worldbank.org/index.php/catalog/1970},
#' accessed 5 November 2020). The data set is ported to R and updated with
#' approximate survey design information derived from the data itself.
#' The units are households and the variables are expenditures on main categories,
#' poverty measures and structural information including weights and sample design.
#'
#' Absolute poverty measures use a poverty line of Lek 4891 (2002 prices).
#' Extreme poverty measures use a poverty line where the basic nutritional needs are
#' difficult to meet.
#' The headcount poverty variable is an indicator for the income of the household \eqn{y_i}
#' being below the (absolute or extreme) poverty line \eqn{z}.
#' The poverty gap variable measures the relative distance to the poverty line: \eqn{(z-y_i)/z}.
#' The poverty depth variable is the square of the poverty gap variable, i.e. \eqn{[(z-y_i)/z]^2},
#' giving more weight to the poorer among the poor and thus describing the inequality
#' among the poor.
#'
#' The survey design is a stratified clustered two stage design.
#' The primary sampling units are enumeration zones.
#' The strata are the crossing of prefecture and urbanicity and the allocation of the
#' psu sample to the strata is proportional to the number of households.
#' Within strata the psu are sampled with probability proportional to number of households.
#' Within psu a simple random sample of 8 households was selected.
#' The weights are calibrated to population margins.
#' All survey design informations except the strata and the weights are approximated
#' through the weights using assumptions on the design.
#' Since the data set has undergone data protection measures and the survey design
#' is approximate only, inference to the population does not yield exact results.
#' However, the complexity of the data and of the survey design are realistic.
#'
#' The size of the household is not on the original data set.
#' However, the transformation \code{capita <- round(0.07527689 * totcons/rcons, 0)}
#' yields the number of persons in the household.
#'
#' @format A data frame with 6671 rows and 26 variables
#' \describe{
#'    \item{ psu }{ primary sampling unit (psu) }
#'    \item{ hhid }{ unique household identifier (100*psu+hh) }
#'    \item{ hh }{ household number per psu }
#'    \item{ prefectu }{ prefecture }
#'    \item{ urban }{ urbanicity (Urban=1, Rural=2) }
#'    \item{ strat }{ stratum }
#'    \item{ region }{ region }
#'    \item{ totcons }{ total consumption of hh }
#'    \item{ rcons }{ real mean per capita consumption }
#'    \item{ rfood }{ real food consumption per capita }
#'    \item{ rtotnfoo }{ real non food consumption per capita }
#'    \item{ reduexpp }{ real education consumption per capita }
#'    \item{ rdurcons }{ real durable consumption per capita }
#'    \item{ rtotutil }{ real utilities consumption per capita}
#'    \item{ egap0 }{ extreme headcount poverty }
#'    \item{ egap1 }{ extreme poverty gap }
#'    \item{ egap2 }{ extreme poverty depth }
#'    \item{ agap0 }{ absolute headcount poverty }
#'    \item{ agap1 }{ absolute poverty gap }
#'    \item{ agap2 }{ absolute poverty depth }
#'    \item{ weight }{ final cross-sectional weight }
#'    \item{ nph }{ number of psu in stratum population }
#'    \item{ mph }{ number of households in stratum population }
#'    \item{ mphi }{ number of households in sampled psu }
#'    \item{ pi1 }{ psu inclusion probability }
#'    \item{ pi2 }{ household inclusion probability }
#'  }
#'
#' @references
#' \url{http://www.instat.gov.al/media/1547/lsms_2012_eng.rar}
#'
#' @note
#' With R package \code{\link{survey}} a survey design object can be built with, e.g., \code{svydesign(~psu + hhid , strata= ~strat, fpc= ~pi1 +pi2,  weight= ~weight, data=lival, pps="brewer")}.
#'
#' @examples
#' library(survey)
#' data(lival)
#' lival$capita <- with(lival, round(0.07527689 * totcons / rcons, 0))
#' lival.des <- svydesign(~psu + hhid , strata= ~strat, fpc= ~pi1 +pi2,
#'                       weight= ~weight, data=lival, pps="brewer")
#' svymean(~totcons, lival.des, deff=TRUE)
"lival"
