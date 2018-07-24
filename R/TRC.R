TRC <-
function(data, weights, overlap=3, mincor=0, 
         robust.regression="rank", gamma=0.5, 
         prob.quantile=0.75, alpha=0.05, md.type="m",
         monitor=FALSE)
{
# Multivariate Outlier Detection in Survey Data
# TRC algorithm as described in:
# B\'eguin, C. and Hulliger B., (2002),
# EUREDIT Workpackage x.2 D4-5.2.1-2.C
# Develop and evaluate new methods for statistical outlier 
# detection and outlier robust multivariate imputation,
# Technical report, EUREDIT 2002.
# 
# Program by C. B\'eguin
# Modification : 27 March 2003 (Beat Hulliger) (R-Version of TRC030313.ssc)
# Modification : 22 December 2005 (Beat Hulliger) Correct standardization of rank correlations if no imputation
# Modification : 16 February 2006 (Beat Hulliger) Correct ordering if p=2
# Modification : 30 March 2006 (Beat Hulliger) default md.type="m", default robust.regression="rank"
# 22.8.2014: Output as result of function
#  Note; The data is not standardised to unit scales.
# Copyright : Swiss Federal Statistical Office, 2002
# gamma: minimal overlap of variables for regression
# robust.regression: type of regression "irls" or based on rank correlation
# alpha: Quantile for F-distribution for cut-off
# md.type=="c" conditional, md.type=="m" marginal Mahalanobis-distance
#
############ Computation time start ############
#
	calc.time <- proc.time()	
#
############ Preprocessing ############
#
	if(!is.matrix(data)) data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	if (alpha<0.5) alpha <- 1-alpha
	if (missing(weights)) (weights <- rep(1,n))
	if (length(weights)!=n) stop("\n Sampling weights of wrong length")
	new.indices <- which(apply(is.na(data),1,prod)==0)
	discarded <- NA
	nfull <- n
	if (length(new.indices)<n) 
	{
		discarded <- which(apply(is.na(data),1,prod)==1)
		cat("Warning: missing observations",discarded,"removed from the data\n")
		data  <-  data[new.indices,]
		weights <- weights[new.indices]
		n <- nrow(data)
	}
	missing.matrix <- 1-is.na(data)
	cat("\n Number of missing items: ", sum(is.na(data)), ", percentage of missing items: ", mean(is.na(data))," \n")
	if (monitor) cat("End of preprocessing in",proc.time()-calc.time,"seconds \n")
#
############ Estimation of location and scatter ############
#
	if (monitor) spearman.time <- proc.time()		
	medians <- apply(data, 2, weighted.quantile, w = weights)
	correction.constant.mad <- 1/qnorm(0.75)
	mads <- correction.constant.mad*apply(abs(sweep(data, 2, medians, "-")), 2, weighted.quantile, w = weights)
	if(sum(mads <= 0) > 0) 
	{
		cat("Some mads are 0. Using", prob.quantile, "quantile absolute deviations!\n")
		mads <- (1/qnorm(0.5*(1+prob.quantile)))*apply(abs(sweep(data, 2, medians, "-")), 2, weighted.quantile, w = weights, prob = prob.quantile)
		if(sum(mads <= 0) > 0)
		{
			cat("The following variable(s) have", prob.quantile, "quantile absolute deviations equal to 0 :",which(mads == 0),"\n")
			stop("Remove these variables or increase the quantile probablity\n")
		}
	}
	if (prod(missing.matrix)==1) # No missing values
	{		
		weighted.ranks <- matrix(0,n,p)
		for (i in 1:p)
		{
			weighted.ranks[,i] <- (apply(data[,i,drop=FALSE],1,.sum.weights,weights=weights,observations=data[,i])
										+0.5*apply(data[,i,drop=FALSE],1,.sum.weights,weights=weights,observations=data[,i],lt=FALSE)
										 +0.5)
		}
		scatter <- (12*(t(weighted.ranks)%*%(weights*weighted.ranks)) / (t(missing.matrix)%*%(weights*missing.matrix))^3-3)
		scatter[scatter>1] <- 1
		scatter[scatter<(-1)] <- -1
            scatter <- 2*sin(pi*scatter/6)
		if (monitor) 
		{
			cat("Spearman Rank Correlations (truncated and standardized):\n")
			print(scatter)
			cat("End of Spearman rank correlations estimations in",proc.time()-spearman.time,"seconds\n")
		}
		if (monitor) 
		{
			cat("No imputation\n")
		}
	}
	else
	{
		scatter <- matrix(0,p,p)
		size.of.cor.sets <- t(missing.matrix)%*%missing.matrix
		if (sum(size.of.cor.sets<overlap)>0)
		{
			cat("Warning: ",(sum(size.of.cor.sets<overlap)-p)/2," couples of variables have less than ",overlap,
			 		" observations in common, therefore their rank correlations will be set to 0.\n")
		}
		if (monitor) cat("Computing Spearman Rank Correlations :\n")
		for (i in 1:(p-1))
		{
			if (monitor) cat("i=",i,"\n")
			for (j in (i+1):p)
			{
				if (monitor) cat(" j=",j,"\n")
				if (size.of.cor.sets[i,j]>=overlap)
				{
					common.observations <- missing.matrix[,i]&missing.matrix[,j]
					weighted.ranks.i <- (apply(data[common.observations,i,drop=FALSE],1,.sum.weights,weights=weights[common.observations],observations=data[common.observations,i])
												+0.5*apply(data[common.observations,i,drop=FALSE],1,.sum.weights,weights=weights[common.observations],observations=data[common.observations,i],lt=FALSE)
										 			+0.5)
					weighted.ranks.j <- (apply(data[common.observations,j,drop=FALSE],1,.sum.weights,weights=weights[common.observations],observations=data[common.observations,j])
												+0.5*apply(data[common.observations,j,drop=FALSE],1,.sum.weights,weights=weights[common.observations],observations=data[common.observations,j],lt=FALSE)
										 			+0.5)
					scatter[i,j] <- 12*sum(weights[common.observations]*weighted.ranks.i*weighted.ranks.j)/sum(weights[common.observations])^3-3
				}
				
			}		}
		scatter <- scatter + t(scatter) + diag(p)
		scatter[scatter>1] <- 1
		scatter[scatter<(-1)] <- -1
		scatter <- 2 * sin(pi * scatter / 6) # Standardization put before imputation 13.03.03 Beat Hulliger
		if (monitor) 
		{
			cat("Spearman Rank Correlations (truncated and standardized):\n")
			print(scatter)
			cat("End of Spearman rank correlations estimations in",proc.time()-spearman.time,"seconds\n")
		}
#
######## Ad hoc imputation of missing values ##############
#
		imputation.time <- proc.time()
		variables.to.be.imputed <- which(apply(missing.matrix,2,prod)==0)
		# regressors correlations with small support are set to 0
    regressors.cor <- (scatter-diag(p))[variables.to.be.imputed,]*(size.of.cor.sets[variables.to.be.imputed,]>=(gamma*n)) 
		regressors.cor <- as.matrix(t(regressors.cor))
		if (monitor) cat("Regressors correlations\n", regressors.cor)
		if (length(variables.to.be.imputed)>1) regressors.list.ordered <- apply(-abs(regressors.cor),1,order) else 
                  regressors.list.ordered <- as.matrix(t(order(-abs(regressors.cor))))
		for (v in 1:length(variables.to.be.imputed))
		{
			observations.to.be.imputed  <- (!missing.matrix[,variables.to.be.imputed[v]])
			if (monitor) cat("Variable",variables.to.be.imputed[v],":\n")
      # loop over candidate regressors
      r <- 0
			repeat
			{
				r <- r+1
				if (abs(regressors.cor[v,regressors.list.ordered[r,v]])<mincor)
				{
          # First in regressors.list.ordered is largest cor
					cat("No eligible regressor found for variable",v,"observation(s)",which(observations.to.be.imputed),".\n Try to relax the regressor eligibility conditions.\n")
					stop()
				}
				if (sum(observations.to.be.imputed & missing.matrix[,regressors.list.ordered[r,v]])==0) next
				observations.imputed.by.r.on.v <- which(observations.to.be.imputed & missing.matrix[,regressors.list.ordered[r,v]])
				k <- length(observations.imputed.by.r.on.v)
				common.observations <- missing.matrix[,variables.to.be.imputed[v]]&missing.matrix[,regressors.list.ordered[r,v]]
				if (robust.regression=="irls") {
					regression.coeff <- rlm(data[common.observations,variables.to.be.imputed[v]]~data[common.observations,regressors.list.ordered[r,v]],weights=weights[common.observations])$coefficients
				} else {
					regression.coeff <- c(0,regressors.cor[v,regressors.list.ordered[r,v]]*mads[variables.to.be.imputed[v]]/mads[regressors.list.ordered[r,v]])
					regression.coeff[1] <- medians[variables.to.be.imputed[v]]-regression.coeff[2]*medians[regressors.list.ordered[r,v]]
				}
				data[observations.imputed.by.r.on.v,variables.to.be.imputed[v]] <- 
                                    matrix(c(rep(1,k),data[observations.imputed.by.r.on.v,regressors.list.ordered[r,v]]),k,2) %*%regression.coeff
				observations.to.be.imputed[observations.imputed.by.r.on.v] <- FALSE
				if (monitor) cat(" ",k,"observations imputed using regressor",regressors.list.ordered[r,v] ,
					            "(cor=",scatter[variables.to.be.imputed[v],regressors.list.ordered[r,v]],
					            "slope=",regression.coeff[2],
					            "intercept=",regression.coeff[1],")\n")
				if (sum(observations.to.be.imputed)>0) next
				break		
			}			
		}
		if (monitor) cat("End of imputation in",proc.time()-imputation.time,"seconds\n")			
	}
	#scatter <- 2*sin(pi*scatter/6) #standardization moved in front of imputation
	scatter <- t(t(mads * scatter) * mads) 
	new.basis <- eigen(scatter)$vectors
	data <- data %*% new.basis 
	center <- apply(data, 2, weighted.quantile, w = weights)
	scatter <- (correction.constant.mad*apply(abs(sweep(data, 2, center, "-")), 2, weighted.quantile, w = weights))^2
	if(sum(scatter == 0) > 0) 
	{
		cat("Some mads are 0. Using", prob.quantile, "quantile absolute deviations!\n")
		scatter <- ((1/qnorm(0.5*(1+prob.quantile)))*apply(abs(sweep(data, 2, center, "-")), 2, weighted.quantile, w = weights, prob = prob.quantile))^2
		if(sum(scatter == 0) > 0)
		{
			stop("Please, increase the quantile probability\n")
		}
	}
	center <- as.vector(new.basis %*% center)
	scatter <- as.matrix(new.basis %*% diag(scatter) %*% t(new.basis))
	data <- data %*% t(new.basis)
#
############ Mahalanobis distances ############
#  
#	dist.with.imputed.values <- mahalanobis(data, center, scatter)
#              var.with.imputed.values <- var(data)
	data[!missing.matrix] <- NA
	s.patterns <- apply(matrix(as.integer(is.na(data)),n,p),1,paste,sep="",collapse="")
	perm <- order(s.patterns)
	data <- data[perm,]
	s.patterns <- s.patterns[perm]
	s.counts <- as.vector(table(s.patterns))
	s.id <- cumsum(s.counts)
	S <- length(s.id)
	missing.items <- is.na(data[s.id,,drop=FALSE])
	nb.missing.items <- apply(missing.items,1,sum)
	indices <- (!missing.items[1,])
	if (md.type=="c") metric <- solve(scatter) else metric <- scatter
	dist <- mahalanobis(data[1:s.id[1],indices,drop=FALSE],center[indices],metric[indices,indices],inverted=(md.type=="c"))*p/(p-nb.missing.items[1])
		if (S>1)
		{
			for (i in 2:S)
			{
				indices <- (!missing.items[i,])
				dist <- c(dist,mahalanobis(data[(s.id[i-1]+1):s.id[i],indices,drop=FALSE],center[indices],
				                             metric[indices,indices,drop=FALSE],inverted=(md.type=="c"))*p/(p-nb.missing.items[i]))
			}
		}
#
############ Choice of outliers ############
#  
# Nominate the outliers using the original numbering (without discarded obs.)
#
	cutpoint <- qf(alpha,p,n-p)/qf(0.5,p,n-p)*median(dist)
	dist <- dist[order(perm)]
	good <- (1:n)[dist < cutpoint] 
	outliers <- (1:n)[ - good]	
# outliers and distances in original numbering with full dataset
      outn <- logical(n)
	  outn[outliers] <- TRUE
	  outnfull <- logical(nfull)
	  outnfull[new.indices] <- outn
	  distnfull <- rep(NA,nfull)
	  distnfull[new.indices] <- dist
#
############ Computation time stop ############
#
	calc.time <- proc.time() - calc.time
#
############ Results ############
#
	TRC.r <- list(sample.size = n, 
				number.of.variables = p, 
				number.of.missing.items=nb.missing.items,
				significance.level = alpha, 
				computation.time = calc.time, 
				medians=medians, mads=mads,
				center = center, 
				scatter = scatter, 
				robust.regression=robust.regression,
				md.type=md.type,
                cutpoint=cutpoint)
#
############ Output ############
#
	cat("\n", "TRC has detected", length(outliers), "outlier(s) in", calc.time[1], "seconds.\n\n")
return(invisible(list(output=TRC.r,outind = outnfull, 
                      dist = distnfull)))
}

