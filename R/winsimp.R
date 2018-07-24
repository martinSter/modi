Winsimp <-
function(data,center,scatter,outind,seed=1000003)
# Imputation under the multivariate normal model after winsorization
# (robust) center and scatter must be given as arguments 
# Beat Hulliger
# 22.5.2009, 7.8.2009
# 22.8.2014 Output as function argument


{
############ Computation time start ############
#
	calc.time <- proc.time()
	
outind<-as.logical(outind)
# Mahalanobis distance (not squared)
data.wins <- as.matrix(data)
MD<-sqrt(MDmiss(data.wins,center,scatter))
cutpoint<-min(MD[outind],na.rm=TRUE)
# robustness weight
u <- ifelse(MD <=cutpoint,1,cutpoint/MD)
#  winsorization for outliers (only outliers!)
data.wins[outind,] <- as.matrix(sweep(sweep(sweep(data[outind,],2,center,'-'),
                                               1,u[outind],'*'),
									2,center,'+'))
# imputation for missing values
rngseed(seed)
s <- prelim.norm(data.wins)
data.imp<-imp.norm(s,makeparam.norm(s,list(center,scatter)),data.wins)

if (sum(is.na(data.imp))>0) cat("There are missing values in the imputed data set.\n")
#
############ Computation time stop ############
#
	calc.time <- proc.time() - calc.time
#
############ Results ############
#

Winsimp.r <-  list(cutpoint=cutpoint, proc.time = calc.time, n.missing.before=sum(is.na(data)),n.missing.after=sum(is.na(data.imp)))
#cat("Results are in winsimp.r and the imputed data is in winsimp.i")

return(invisible(list(output=Winsimp.r,imputed.data=data.imp)))

}

