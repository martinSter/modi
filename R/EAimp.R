EAimp <-
function(data, weights , outind, reach="max",
         transmission.function = "root", power=ncol(data), distance.type = "euclidean",
         duration=5, maxl = 5, kdon=1, monitor = FALSE, 
         threshold=FALSE, deterministic=TRUE, fixedprop=0 )
{
# EPIDEMIC Algorithm for Multivariate Impuation 
#
# B\'eguin, C. and Hulliger, B. (2004) Multivariate outlier detection in incomplete survey data: 
# the epidemic algorithm and transformed rank correlations, JRSS-A, 167, Part 2, pp. 275?294.
#
# Program by Cédric Béguin and Beat Hulliger 
# Created : Wednesday, January 24, 2001
# Modifications : 4 August 2009 Beat Hulliger
# 4.7.2003: Conversion to R from EA030313.ssc by Beat Hulliger
# 27.3.2009: Modular programming and packaging: Beat Hulliger 
# 23.8.2014: Output passed as function results
#  
# Copyright Swiss Federal Statistical Office and EUREDIT 2001-2006, FHNW 2007-2009
# 
# discrete: if TRUE changes the correction for missing values (instead of *(p/q) +(p-q) when summing (corresponds to corr=0.5)
# reach: Transmission.distance. 
# maxl: Maximum number of steps without change (finish)
# threshold: Infect all points with probability above 1-0.5^(1/maxl)
# fixedprop: Fixed proportion to be infected
# deterministic: Infect points with largest prob. (expected number)
# outind: a logical vector with TRUE if an outlier.
#
#
############ Computation time start ############
#
	calc.time <- proc.time()[1]	
############ Dimensions ############
#
# use all data: complete non-responses are to be replaced 
#	completely (similar to outliers) 
#
	n <- nrow(data)
	p <- ncol(data)
	if (missing(weights)) weights<-rep(1,n)
  if (missing(outind)) cat("No outlier indicator is given\n")
  if (sum(is.na(outind))>0) {
    cat("Missing values in outlier indicator set to FALSE.\n")
    outind[is.na(outind)] <- FALSE
  }
	response<-!is.na(data)
	complete.records <- apply(response, 1, prod)
	usable.records <- apply(response, 1, sum) >= p/2
	cat("\n Dimensions (n,p):", n, p)
	cat("\n Number of complete records ", sum(complete.records))
	cat("\n Number of records with maximum p/2 variables missing ", sum(usable.records))
#
# Standardization of weights to sum to sample size
#
	np <- sum(weights)
	weights <- as.single((n * weights)/np)	
#
#
# Calculation of distances (independently of EAdet)
EA.dist.res<-.EA.dist(data, n=n,p=p,weights = weights,reach=reach,
                      transmission.function = transmission.function, power=power, distance.type = distance.type, 
                      maxl = maxl)
if (monitor) cat("\n\n Distances finished")

if (length(EA.dist.res$counterprobs)!=n*(n-1)/2) cat("\n Distances not on same number of observations")
#
# Imputation with reverse epidemic and random imputation among nearest potential donors
#
# Preparation
#
  imp.data<-data
  to.impute<-apply(response,1,prod)==0 | outind
  imp.ind<-which(to.impute)
  n.imp<-length(imp.ind)
  cat("\n Number of imputands is ",n.imp)
  # reach d0 is set to maximum in order to reach all outliers (max.min.di of .EA.dist)
  d0 <- reach
  cat("\n Reach for imputation is ",d0)
#
# Start of imputations
#
  for (i in 1:n.imp) 
  {
# Start epidemic at the observation that needs imputation
    imputand<-imp.ind[i]
    imp.time<-1
    imp.infected<-rep(FALSE,n)
    imp.infected[imputand]<-TRUE
    if (outind[imputand]) missvar<-(1:p) else missvar<-which(!response[imputand,])
    n.missvar<-length(missvar)  
    # donors must have values for the missing ones and must not be an outlier
    potential.donors<-apply(as.matrix(response[,missvar]),1,sum)==n.missvar & !outind
    potential.donors[imputand]<-FALSE
	if (monitor) {
		cat("\n Imputand: ", imputand)
		cat(" missvar: ",missvar)  
    		cat(" #pot. donors: ",sum(potential.donors))
	}
    # if no donor can impute all missing values of the imputand then maximise number of imputed values
    if (sum(potential.donors)==0) 
    {
      max.miss.val<-max(apply(as.matrix(response[-imputand,missvar]),1,sum))
      cat("\n No common donor for all missing values. ", n.missvar-max.miss.val," missing values will remain")
      potential.donors<-apply(response[,missvar],1,sum)==max.miss.val
    }
    donors<-rep(FALSE,n)
    new.imp.infected <- imp.infected
    n.imp.infected <- sum(imp.infected)
    hprod <- rep(1, n)
    imp.infection.time <- rep(0, n)
    imp.infection.time[imp.infected] <- imp.time
#
############ Run epidemic from imputand until at least one potential donor is infected
#  
	repeat {
            if (sum(donors)>=kdon | imp.time>duration) break		
            #cat("\n imp.time = ", imp.time, " , imp.infected = ", n.imp.infected)
		#print(memory.size())
		imp.time <- imp.time + 1
		old.imp.infected <- imp.infected
		if(sum(new.imp.infected) > 1) {
		 hprod[!imp.infected] <- hprod[!imp.infected] *  apply(sweep(sweep(matrix(EA.dist.res$counterprobs[apply(as.matrix(which(!imp.infected)),
				1, .ind.dijs, js = which(new.imp.infected), n = n)], sum(new.imp.infected), n - n.imp.infected),1,
				weights[new.imp.infected],"^" ),2,weights[!imp.infected],"^"), 2, prod)
					
		}
		else {
			if(sum(new.imp.infected) == 1)
				  hprod[!imp.infected] <- hprod[!imp.infected] * EA.dist.res$counterprobs[apply(as.matrix(which(!imp.infected)), 1, 
				  .ind.dijs, js = which(new.imp.infected), n = n)]^(weights[new.imp.infected]* 
					weights[!imp.infected])
			
		}
		if (deterministic) {
			n.to.infect <- max(1,round(sum(1 - hprod[!imp.infected]))) # At least 1 infection at each step
                  # Rank is maximum to allow infection
			imp.infected[!imp.infected] <- rank(1 - hprod[!imp.infected],ties.method="max")>=n-n.imp.infected-n.to.infect
		} else {
			if (threshold) {imp.infected[!imp.infected] <- hprod[!imp.infected]<=0.5^(1/maxl)} else {
				if (fixedprop>0) {
					n.to.infect <- max(1,floor(fixedprop* sum(!imp.infected)))
					imp.infected[!imp.infected] <- rank(1 - hprod[!imp.infected])>=n-n.imp.infected-n.to.infect
   				} else 
		      imp.infected[!imp.infected] <- as.logical(rbinom(n - n.imp.infected, 1, 1 - hprod[!imp.infected]))}}
		new.imp.infected <- imp.infected & (!old.imp.infected)
		n.imp.infected <- sum(imp.infected)
		imp.infection.time[new.imp.infected] <- imp.time
           donors<-potential.donors & imp.infected
 		if (monitor) cat(", donors: ", which(donors))
		next		
	}
    if (imp.time>duration) cat("No donor found") else 
    {
      if (sum(donors)>1) don.imp<-sample(which(donors),1) else don.imp<-which(donors)
      if (monitor) cat(". #donors: ", sum(donors),", chosen donor: ",don.imp)
      imp.data[imputand,missvar]<-data[don.imp,missvar]
    }
    #cat("\n memory use",memory.size())
}
calc.time <- round(proc.time()[1] - calc.time, 5)
cat("\n\n Number of remaining missing values is ", sum(is.na(imp.data)))
#
############ Results ############
#
	EAimp.r <- list(sample.size = n, number.of.variables = p, n.complete.records = sum(complete.records), 
		n.usable.records = sum(usable.records),duration=duration,reach=d0, threshold=threshold, 
		deterministic=deterministic, computation.time = calc.time)
return(invisible(list(output=EAimp.r,imputed.data=imp.data)))
#
############ Output ############
#
	cat("\n", "EA imputation has finished in", calc.time, "seconds.", "\n")
#		cat("The results are in EAimp.r and in EAimp.data \n")
}

