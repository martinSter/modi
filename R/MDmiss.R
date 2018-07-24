MDmiss<-function(data,center,cov)
# Mahalanobis distance for data with missing values
# The center and scatter must be complete
# The function loops over the observations
# Beat Hulliger
# 26.8.2014
{
  if(!is.matrix(data)) data<-as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  missings<-is.na(data)
#  if (sum(apply(missings,1,prod)==1)>0 )
#  {
#    cat("Remove completely missing observations!\n")
#    return(NA)
#  }
  # count number of responded variables per observation
  resp.dim<-apply(!missings,1,sum)
  # center the data
  data<-sweep(data, 2, center)
  mdm<-numeric(n)
  for (i in 1:n) {
    if (resp.dim[i]>0){
    resp<-!missings[i,]
    x<-data[i,resp]
    mdm[i]<-t(x) %*% solve(cov[resp,resp]) %*% x
    } else mdm[i]<-NA
  }
  # correction for number of responding variables
  mdm<-mdm*p/resp.dim
  return(mdm)
}