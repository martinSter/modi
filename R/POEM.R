POEM <-
function(data,weights,outind,errors,missing.matrix,alpha=0.5,beta=0.5,reweight.out=FALSE,c=5,
                 preliminary.mean.imputation=FALSE,monitor=FALSE)
{
# POEM Algorithm for multivariate weighted imPutation for Outliers, Edit failure and Missing values
#
# POEM algorithm as described in:
# B?guin, C. and Hulliger B., (2002),
# EUREDIT Workpackage x.2 D4-5.2.1-2.C
# Develop and evaluate new methods for statistical outlier 
# detection and outlier robust multivariate imputation,
# Technical report, EUREDIT 2002.
#
# Program by C?dric B?guin
# Last modified : 7 August 2009 (Beat Hulliger)
# Adaptation to R 4.7.2003 (Beat Hulliger from POEM030402.ssc)
# Copyright : Swiss Federal Statistical Office, 2002
# alpha: Weight of failing items
# beta: Condition for link to donor
# c: Tuning constant for redefinition of outliers
### Initial tests ###
#
if (!is.matrix(data)) data<-as.matrix(data)
n <- nrow(data)
p <- ncol(data)
if (missing(weights)) {weights <- rep(1,n)}
if (is.logical(outind)) outind<-as.numeric(outind)
if (missing(missing.matrix)){missing.matrix <- (1-is.na(data))}
if (missing(errors)){errors <- matrix(1,nrow=n,ncol=p)}
if (!is.vector(weights)) stop("Weights not in vector form","\n")
if (!is.matrix(missing.matrix)) stop("Missing values not in matrix form","\n")
if (!is.vector(outind)) stop("outind not in vector form","\n")	
if (!is.matrix(errors)) stop("Errors not in matrix form","\n")
if (length(weights)!=n) stop("Wrong length of weights")
if (nrow(missing.matrix)!=n | ncol(missing.matrix)!=p) stop("Missing values matrix do not have same dimensions as data","\n")
if (length(outind)!=n) stop("Wrong length of outind","\n")
if (nrow(errors)!=n | ncol(errors)!=p) stop("Errors matrix do not have same dimensions as data","\n")
if (sum(is.na(weights))>0) stop("Missing values in weights","\n")
if (sum(is.na(missing.matrix))>0) stop("Missing values in missing data matrix","\n")
if (sum(is.na(outind))>0) stop("Missing values in outind","\n")
if (sum(is.na(errors))>0) stop("Missing values in errors matrix","\n")
#
# Completely missing
comp.miss<-apply(is.na(data),1,prod)
cat("\n Number of completely missing observations ",sum(comp.miss),"\n")
#
############ Computation time start ############
#
	calc.time <- proc.time()
#
#
### Set all missing values to zero ###
#
old.data <- data
data[!(missing.matrix)] <- 0
#
### Computation of alpha_ij ###
#
alpha.ij <- missing.matrix*(alpha^(1-errors))
good.values <- apply(weights*alpha.ij,2,sum)
old.number.of.nonoutliers <- sum(1-outind)
old.weighted.sum.of.nonoutliers <- sum(weights*(1-outind))
missing.errors <- missing.matrix*errors
#
### Computation of center ###
#
center <- apply((1-outind)*weights*alpha.ij*data,2,sum)/apply((1-outind)*weights*alpha.ij,2,sum)
#
### Centering of data ###
#
data <- sweep(data,2,center,"-") 
#
### Computation of coordinates variances ###
#
variances <- apply((1-outind)*weights*alpha.ij*data^2,2,sum)
variances <- variances/apply((1-outind)*weights*alpha.ij,2,sum)
if (sum(variances==0)>0)
{
	zero.variances <- which(variances==0)
	cat("Warning: Variable(s)",zero.variances,"has (have) zero variance(s)\n")
	stop("\nRemove these variables or reduce the set of outliers\n")
}
#
### Standardization of data ###
#
data <- sweep(data,2,sqrt(variances),"/")
#
### Computation of covariance matrix ###
#
covariance <- (t(alpha.ij*data)%*%((1-outind)*weights*alpha.ij*data))
if (!preliminary.mean.imputation) 
{ # Formula (12) on page 121 of the Euredit Deliverable
	covariance <- covariance/(t(alpha.ij)%*%((1-outind)*weights*alpha.ij))
	if (determinant(covariance)$sign<0) 
	{ 
	cat("Warning: Covariance matrix not positive definite with original data including missing values\n")
	cat("         Choose option preliminary.mean.imputation=T!\n")
	}
} 
# preliminary mean imputation is Formula (13) on page 122 of the Euredit Deliverable.  
# This is equivalent to setting tilde x to zero. And this in turn 
# is equivalent to setting x to the mean before standardising.
 else covariance <- covariance / sum((1-outind)*weights) 
if (monitor) {cat("Covariance matrix of standardised observations\n")
	print(covariance)}
#
### Reweighting of outliers ###
#
if (reweight.out)
{
	MD <- mahalanobis(alpha.ij*data,rep(0,p),covariance)
	#if (!preliminary.mean.imputation) # This condition waived to be compatible with D4-5.2.1-2.C
	{
		MD <- p^2*MD/apply(alpha.ij,1,sum)^2
		if (min(MD)<0) cat ("Warning: Negative Mahalanobis distances\n")
	}
	outind <- 1-(MD > (c*qchisq(2*pnorm(1)-1,p)))
	cat("New set of",sum(outind),"outliers generated\n") 
}
#
### List of observations to be imputed ###
#
observations.with.errors <- (apply(missing.errors,1,sum)<p)
to.be.imputed <- which(outind|observations.with.errors)
imputed <- rep(0,n)
#
### Start of imputation process ###
#
### List of complete and correct donors ###
#
complete.donors <- apply((1-outind)*missing.errors,1,sum)==p
cat("\n Number of complete and error-free observations: ", sum(complete.donors),"\n")
#
for (observation in to.be.imputed) 
{
	#
	### Indicator of potential donnors for the observation ###
	#
	if (outind[observation]==0 & beta<1)
	{
	   potential.donors <- 
					apply((1-outind)*sweep(missing.errors,2,(missing.errors[observation,]),"*"),1,sum) >= beta*p &
					apply(sweep(1-missing.errors,2,(1-missing.matrix[observation,]),"*"),1,sum)==0 &
					apply(sweep(1-missing.errors,2,(1-errors[observation,]),"*"),1,sum)==0								
	}
	else
	{
		potential.donors <- complete.donors
	}
	potential.donors[observation] <- FALSE # Exclude the observation as its own donor	
	switch(sum(potential.donors)+1,
	{# 1
		cat("\n No donor for observation ",observation,"\n")
		cat("All complete error-free observations used as donors.\n")
		potential.donors<-complete.donors
	}, 
	{# 2
		cat("\n Only one donor for observation ",observation,"\n")
	})
	if (monitor) cat("Number of potential donors ",sum(potential.donors),"\n")
#
### Distances of the observation to potential donors ###
#
	distances.to.donors <- mahalanobis(sweep(data[potential.donors,],2,data[observation,],"-")*
                                           sweep(alpha.ij[potential.donors,],2,
                                            alpha.ij[observation,],"*"),rep(0,p),covariance)
	#if (!preliminary.mean.imputation) # This condition waived to be compatible with D4-5.2.1-2.C
		{
			distances.to.donors <- p^2*distances.to.donors/
                                         apply(sweep(alpha.ij[potential.donors,],2,alpha.ij[observation,],"*"),1,sum)^2
		}		
	# 
	### Selection of the donor ###
	#
	min.dist.to.donors <- min(distances.to.donors)
	# donors are the indices for a vector with the indices of the potential donors
	# potential.donors is an indicator vector of length n
	#
	if (is.na(min.dist.to.donors)) donors<-1:sum(potential.donors) else 
		if (min.dist.to.donors<0) stop("Warning: Minimal distance to nearest neighbour negative\n") else
	       donors <- which(distances.to.donors==min.dist.to.donors)
    if (length(donors)==1)
	{
		imputed[observation] <- which(potential.donors)[donors]
	}
	else
	{
		imputed[observation] <- which(potential.donors)[sample(donors,1)]
	}
	if (monitor) {
		cat("Observation",observation,"imputed by donor ",imputed[observation],"\n")
		cat("distance to donor: ",min.dist.to.donors,"\n")
	}
		
}
#
### Imputation ###
#
new.data <- old.data
new.data[as.logical(outind),] <- old.data[imputed[as.logical(outind)],]
non.outliers.errors <- which((1-outind) & observations.with.errors)
new.data[non.outliers.errors,][missing.errors[non.outliers.errors,]<1]  <- 
               old.data[imputed[non.outliers.errors],][missing.errors[non.outliers.errors,]<1]
new.center <- apply(weights*new.data,2,sum)/sum(weights)
new.variances <- apply(weights*sweep(new.data,2,new.center,"-")^2,2,sum)/sum(weights)
#
############ Computation time stop ############
#
	calc.time <- proc.time() - calc.time
#
################## Results ##################
#
   POEM.r <- list(
    preliminary.mean.imputation=preliminary.mean.imputation,
		completely.missing=sum(comp.miss),
    good.values=good.values,
		nonoutliers.before=old.number.of.nonoutliers,
		weighted.nonoutliers.before=old.weighted.sum.of.nonoutliers,
		number.of.nonoutliers.after.reweighting=sum(1-outind),
		weighted.number.of.nonoutliers.after.reweighting=sum(weights*(1-outind)),
		old.center=center, old.variances=variances,
    new.center=new.center, new.variances=new.variances,
    covariance=covariance,
		imputed.observations = to.be.imputed,
		donors = imputed[to.be.imputed],
		outind = outind)
#
################## Output ##################
#
   cat("\n","POEM has imputed",length(to.be.imputed),"observations(s) in",calc.time,"seconds.","\n","\n")
#	cat(" The results are in POEM.r and POEM.i \n")
return(invisible(list(output=POEM.r,imputed.data=new.data)))
}

