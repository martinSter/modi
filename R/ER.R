ER <-
function(data,weights,alpha=0.01,psi.par=c(2,1.25),
         em.steps=100,steps.output=FALSE,Estep.output=FALSE,tolerance=1e-6)
{
##################  Preprocessing of the data  ##################
#
# Removing the unit(s) with all items missing
#
if (!is.matrix(data)) data<-as.matrix(data)
n <- nrow(data)
p <- ncol(data)
if (alpha<0.5) alpha<-1-alpha
if (missing(weights)) weights <- rep(1,n)
new.indices <- which(apply(is.na(data),1,prod)==0)
if (length(new.indices)<n) 
{
    cat("Warning: missing observations",which(apply(is.na(data),1,prod)==1),"removed from the data\n")
    data <- data[new.indices,]
    weights <- weights[new.indices]
    n <- nrow(data)
}
if (steps.output) cat("End of preprocessing\n")
#
############ Computation time start ############
#
    calc.time <- proc.time()
#
# Order the data by missingness patterns : 
# s.patterns = vector of length n, with the missingness patterns stocked as strings of the type "11010...011" with "1" for missing.
# data, weights and s.patterns ordered using s.patterns' order
#
s.patterns <- apply(matrix(as.integer(is.na(data)),n,p),1,paste,sep="",collapse="")
perm <- order(s.patterns)
data <- data[perm,]
s.patterns <- s.patterns[perm]
weights <- weights[perm]
#
# Missingness patterns stats :
#
# s.counts = counts of the different missingness patterns ordered alphabetically.
# s.id = indices of the last observation of each missingness pattern in the dataset ordered by missingness pattern.
# S = total number of different missingness patterns
# missing.items = missing items for each pattern
# nb.missing.items = number of missing items for each pattern
#
s.counts <- as.vector(table(s.patterns))
s.id <- cumsum(s.counts)
S <- length(s.id)
missing.items <- is.na(data[s.id,,drop=FALSE])
nb.missing.items <- apply(missing.items,1,sum)
if (steps.output) cat("End of missingness statistics\n")
#
#
################## Preparation for call to ER ##################
#
# This is the old BEM step 1 for the choice of the initial good subset
# 
#
        #
        # mean and covariance matrix computed by ER 
        #
            if (nb.missing.items[1]==0)
            {
            #
            # Case where some observations are complete (no missing item) => 
            # start.mean and start.var are calculated on the complete observations
            #
                cov.complete <- cov.wt(data[1:s.counts[1],],wt=weights[1:s.counts[1]]/sum(weights[1:s.counts[1]]))
                mean.start <- cov.complete$center
                var.start <- cov.complete$cov
            }
            else
            {
            #
            # Case where all observations have missing items
            #
               mean.start <- apply(data,2,weighted.mean,w=weights,na.rm=TRUE)
               var.start <- diag(apply(data,2,weighted.var,w=weights,na.rm=TRUE))
            }
		if (steps.output) {
      cat("\n","start.mean: ",mean.start,"\n","start.var: ")
      print(var.start)
      cat("\n")
		}
            ER.result <- .ER.normal(data=data,weights,psi.par=psi.par,np=sum(weights),p=p,s.counts=s.counts,s.id=s.id,S,
                                    missing.items=missing.items, nb.missing.items,
                                    start.mean=mean.start,
                                    start.var=var.start,
                                    numb.it=em.steps,Estep.output=Estep.output,
                                    tolerance=tolerance)

            ER.mean <- ER.result$theta[1,2:(p+1)]
            ER.var <- ER.result$theta[2:(p+1),2:(p+1)]
#
############ Computation time stop ############
#
    calc.time <- proc.time() - calc.time
#
################## Step 5 ##################
#
# Nominate the outliers using the original numbering
#
    dist <- ER.result$dist 
    good <- dist <= qchisq(alpha,p)
    outliers <- perm[!good]
#
#
################## Output ##################
#
   cat("\n","ER has detected",sum(!good),"outlier(s) in",calc.time,"seconds.","\n","\n")
   if (!ER.result$convergence) cat("\n","ER did not converge.","\n")
################## Results ##################
#
return(list(sample.size = n,
            number.of.variables = p,
            significance.level = alpha,
            computation.time = calc.time,
            good.data = perm[good],
            outliers = outliers,
            center = ER.mean,
            scatter = ER.var,
            dist = dist[order(perm)],
            rob.weights = ER.result$rob.weights))  
}

