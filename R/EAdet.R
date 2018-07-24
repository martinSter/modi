EAdet <-
function(data, weights ,reach="max", transmission.function = "root", power=ncol(data), distance.type = "euclidean", 
	maxl = 5, plotting = TRUE, monitor = FALSE, prob.quantile = 0.9, 
	random.start = FALSE, fix.start, threshold=FALSE, deterministic=TRUE, rm.missobs=FALSE, verbose=FALSE)
{
# EPIDEMIC Algorithm for Multivariate Outlier Detection
#
# B\'eguin, C. and Hulliger, B. (2004) Multivariate outlier detection in incomplete survey data: 
# the epidemic algorithm and transformed rank correlations, JRSS-A, 167, Part 2, pp. 275?294.
#
# Program by C\'edric B\'eguin and Beat Hulliger 
# Created : Wednesday, January 24, 2001
# Last modification : 4 August 2009 Beat Hulliger
# Conversion to R from EA030313.ssc by Beat Hulliger (4.7.2003)
# Modular programming and packaging: Beat Hulliger 27.3.2009
# 22.8.2014 Output as objects
# Copyright Swiss Federal Statistical Office and EUREDIT 2001-2006, FHNW 2007-2009
############ Dimensions ############
#
#

#
if (!is.matrix(data)) data<-as.matrix(data)
n <- nrow(data)
p <- ncol(data)
if (missing(weights)) weights <- rep(1,n)
# Removing the unit(s) with all items missing
new.indices <- which(apply(is.na(data),1,prod)==0)
miss.indices<-apply(is.na(data),1,prod)==1
missobs<-which(miss.indices)
discarded<-NA
nfull<-n
if ((length(new.indices)<n) & rm.missobs)
{
	discarded<-missobs
    cat("Warning: missing observations",discarded,"removed from the data\n")
    data <- data[new.indices,]
    weights <- weights[new.indices]
    n <- nrow(data)
}
	complete.records <- apply(!is.na(data), 1, prod)
	usable.records <- apply(!is.na(data), 1, sum) >= p/2
	if (verbose) {cat("\n Dimensions (n,p):", n, p)
	cat("\n Number of complete records ", sum(complete.records))
	cat("\n Number of records with maximum p/2 variables missing ", sum(usable.records),"\n")
	}
  power <- as.single(power)
#
# Standardization of weights
#
	np <- sum(weights)
	weights <- as.single((n * weights)/np)	#
#
############ Computation time start ############
#
	calc.time <- proc.time()[1]	#
############ Calibraton and setup ############
#
	medians <- apply(data, 2, weighted.quantile, w = weights, prob = 0.5)
	data <- sweep(data, 2, medians, "-")
	mads <- apply(abs(data), 2, weighted.quantile, w = weights, prob = 0.5)
	qads <- apply(abs(data), 2, weighted.quantile, w = weights, prob = prob.quantile)
	if(sum(mads == 0) > 0) {
		cat("\n Some mads are 0. Standardizing with ", prob.quantile, " quantile absolute deviations!")
		
		if(sum(qads == 0) > 0)
			cat("\n Some quantile absolute deviations are 0. No standardization!")
		else data <- sweep(data, 2, qads, "/")
	}
	else data <- sweep(data, 2, mads, "/")
#
# Calculation of distances
	EA.dist.res<-.EA.dist(data, n=n,p=p,weights = weights,reach=reach,
				transmission.function = transmission.function, power=power, distance.type = distance.type, 
				maxl = maxl)
	if (monitor) cat("\n\n Distances finished")
# The distances calculated by EA.dist are the counterprobabilities in single precision. 
if (monitor) {
cat("\n memory use after distances: ",memory.size())
cat("\n Index of sample spatial median is ",EA.dist.res$output[1])
cat("\n Maximal distance to nearest neighbor is ", EA.dist.res$output[2])
cat("\n Transmission distance is ", EA.dist.res$output[3], "\n")
}
#
############ Initialisation ############
# 
	if (verbose) cat("\n\n Initialisation of epidemic")
	comp.time.init <- proc.time()[1] - calc.time
	if(monitor)
		cat("\n Initialisation time is ", comp.time.init)
	if(random.start)
		start.point <- sample(1:n, 1, prob = weights)
	else {
		if(!missing(fix.start))
			start.point <- fix.start
		else start.point <- EA.dist.res$output[1]
	}
	time <- 1
	infected <- rep(FALSE, n)
	infected[c(start.point)] <- TRUE
	new.infected <- infected
	n.infected <- sum(infected)
	hprod <- rep(1, n)
	
	infection.time <- rep(0, n)
	infection.time[c(start.point)] <- time	#
############ Main loop ############
#  
	repeat {
		if (monitor) cat("\n time = ", time, " , infected = ", n.infected)
		#print(memory.size())
		time <- time + 1
		old.infected <- infected
		if(sum(new.infected) > 1) {
		 hprod[!infected] <- hprod[!infected] *  apply(sweep(sweep(matrix(EA.dist.res$counterprobs[apply(as.matrix(which(!infected)),
				1, .ind.dijs, js = which(new.infected), n = n)], sum(new.infected), n - n.infected),1,
				weights[new.infected],"^" ),2,weights[!infected],"^"), 2, prod)
					
		}
		else {
			if(sum(new.infected) == 1)
				  hprod[!infected] <- hprod[!infected] * EA.dist.res$counterprobs[apply(as.matrix(which(!infected)), 1, 
				  .ind.dijs, js = which(new.infected), n = n)]^(weights[new.infected]* 
					weights[!infected])
			
		}
		if (deterministic) {
			n.to.infect <- sum(1 - hprod[!infected]) # HRK: expected number of infections
			# Do maxl trials for very small inf. prob.
			if (n.to.infect<0.5) n.to.infect <- sum(1 - hprod[!infected]^maxl) 
			n.to.infect<-round(n.to.infect)
			infected[!infected] <- rank(1 - hprod[!infected])>=n-n.infected-n.to.infect
		} else {
			if (threshold) {infected[!infected] <- hprod[!infected]<=0.5^(1/maxl)} else
		    infected[!infected] <- as.logical(rbinom(n - n.infected, 1, 1 - hprod[!infected]))}
		new.infected <- infected & (!old.infected)
		n.infected <- sum(infected)
		infection.time[new.infected] <- time
		if(n.infected == n) {
			break
		}
		if((time - max(infection.time)) > maxl) {
			break
		}
		next		
	}
	duration <- max(infection.time)
#	if(monitor) {
#		last.infection.prob <- 1 - hprod
#	}
if (verbose) 	cat("\n memory use after epidemic: ",memory.size())
# 

############ Computation time stop ############
#
	calc.time <- round(proc.time()[1] - calc.time, 5)

# Default cutpoint
  med.infection.time <- weighted.quantile(infection.time,weights,0.5)
	mad.infection.time <- weighted.quantile(abs(infection.time-med.infection.time),weights,0.5)
  if (verbose) cat("\n med and mad of infection times: ",med.infection.time," and ",mad.infection.time)
	if (mad.infection.time==0) mad.infection.time <- med.infection.time
  cutpoint <- min(med.infection.time+3*mad.infection.time,duration)
  if (verbose) cat("\n Proposed cutpoint is ",min(cutpoint,duration))
# Blowing up to full length
infectedn<-logical(n)
infectedn[infected]<-TRUE
infectednfull<-logical(nfull)
if (nfull>n) infectednfull[new.indices]<-infectedn	else infectednfull<-infectedn	
inf.time<-rep(NA,nfull)
if (nfull>n) inf.time[new.indices]<-infection.time else inf.time<-infection.time

# Set infection time of completely missing to NA
inf.time[!infectednfull]<-NA
# outliers full sample
outlier<-(inf.time>=cutpoint) 
# outlier[is.na(outlier)]<-FALSE
outlier.ind<-which(outlier)
#

####
# Plotting
# not infected are set to high inf.time to show better the healthy on a graph of infection times
if(plotting) 
{
  plot.time<-inf.time
  plot.time[!infectednfull] <- ceiling(1.2 * duration)
  ord <- order(plot.time)
  plot(plot.time[ord],cumsum(weights[ord]), xlab="infection time",  ylab = "(weighted) cdf of infection time")
  abline(v=cutpoint)
}  
#
#
############ Results ############
#
	EAdet.r <- list(sample.size = n, discarded.observations=discarded, missing.observations=missobs,
	             number.of.variables = p, n.complete.records = sum(complete.records), 
		n.usable.records = sum(usable.records), medians = medians, mads = mads, prob.quantile = prob.quantile, 
		quantile.deviations = qads, start = start.point, transmission.function = transmission.function, power=power,
		maxl=maxl,
		max.min.di=EA.dist.res$output[2],	transmission.distance = EA.dist.res$output[3], threshold=threshold, distance.type = distance.type, 
		deterministic=deterministic, number.infected = n.infected, 
                cutpoint=cutpoint, number.outliers=sum(outlier), outliers=outlier.ind,
		duration = duration, computation.time = calc.time, initialisation.computation.time = comp.time.init)


############ Output ############
#
	cat("\n", "EA detection has finished with", n.infected, "infected points in", calc.time[1], "seconds.")
	#	cat("\n The results are in EAdet.r and EAdet.i", "\n")
return(invisible(list(output=EAdet.r,
            infected = infectednfull, 
            infection.time = inf.time, 
            outind=outlier)))
}

