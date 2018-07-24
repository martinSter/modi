GIMCD <- function(data, alpha = 0.05, seedem, seedmcd) {

  # runs em and then mcd
  # Beat Hulliger, 2007
  # 30.10.2015: control of seed for mcd
  ############ Computation time start ############

  calc.time <- proc.time()

  if (missing(seedem)) seedem <- 23456789
  rngseed(seedem)
  if (!is.matrix(data)) data<-as.matrix(data)
  if (alpha<0.5) alpha<-1-alpha
  s <- prelim.norm(data)
  thetahat <- norm::em.norm(s, showits=FALSE)
  imp.data <- norm::imp.norm(s, thetahat, data)

  # MCD algorithm
  if (!missing(seedmcd)) set.seed(seedmcd)
  MCD.cov <- MASS::cov.mcd(imp.data)
  dist <-mahalanobis(imp.data, MCD.cov$center,  MCD.cov$cov)
  n<-nrow(data)
  p<-ncol(data)
  cutpoint <- qf(alpha, p, n-p)*median(dist)/qf(0.5, p, n-p)
  outind<-(dist > cutpoint)

  ############ Computation time stop ############
  #
  calc.time <- proc.time() - calc.time
  #
  ############ Results ############
  #

  cat("GIMCD has detected", sum(outind), "outliers in", calc.time, "seconds.")
  return(list(center=MCD.cov$center, scatter=MCD.cov$cov,
              alpha=1-alpha, computation.time = calc.time[1],
              cutpoint = cutpoint,
              outind=outind, dist=dist))

}
