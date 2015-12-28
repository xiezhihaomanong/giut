giut <- function(..., threshold=0.5) 
# This brings everything together. First, it estimates the proportion of genes
# in each combination of null/alternative hypotheses. As this estimate isn't
# reliable or unbiased, the function then attempts to maximize the p-value by
# searching the parameter space around the estimate. This assumes that the
# estimate is reasonably close, such that the local maximum is conservative.
#
# written by Aaron Lun
# 30 March 2014
{
	all.p <- list(...)
	out <- processPvals(all.p)
	o <- order(out$p.max, decreasing=TRUE) - 1L

	# Initializing the proportions.
	theta0 <- estimateProp(all.p, out$design, threshold=threshold)

	# Computing pmax for each gene. 
	pval <- .Call("R_load_stuff", out$p.max, out$m, theta0, out$design, o, PACKAGE="giut")
	if (is.character(pval)) { stop(pval) }
	return(pval)
}

