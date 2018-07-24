weighted.quantile <-
function(x, w, prob = 0.5,plot=FALSE)
{
#### Weighted quantile function (default weighted median)
#### Automatically removes missing values
#### May plot the empirical distribution function
## Programme by C?dric B?guin and Beat Hulliger
## Copyright 2003 Swiss Federal Statistical Office
		if(missing(w))
			return(quantile(x, prob, na.rm = TRUE))
		else w <- w[!is.na(x)]
		x <- x[!is.na(x)]
		ord <- order(x)
		w <- w[ord]
		x <- x[ord]
		w.ord <- cumsum(w)/sum(w)
		index <- 1:length(x)
		if (plot) plot(x,w.ord,type="s")
		if(min(w.ord)>prob) {
			cat("\n Dominance of one observation!\n")
			lower.k.quant <- 1
		} else lower.k.quant <- max(index[w.ord <= prob])
		upper.k.quant <- min(index[w.ord > prob])
		if(w.ord[lower.k.quant] < prob)
			return(x[upper.k.quant])
		else return((w[lower.k.quant] * x[lower.k.quant] + w[upper.k.quant] * x[upper.k.quant])/(w[
				lower.k.quant] + w[upper.k.quant]))
}

