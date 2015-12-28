estimateProp <- function(pvals, design, threshold=0.5) 
# A function to estimate the proportion of genes for each null/alt combination.
# This assumes that, at the specified threshold, there are no type II errors.
# It then calculates the number of genes with maximum p-values below the threshold
# for each combination of experiments.
{
	lhs <- rhs <- NULL
	for (i in 1:ncol(design)) { 
		curnull <- design[,i]==1L
		vested <- colSums(design[curnull,,drop=FALSE])
		p.max <- do.call(pmax, pvals[curnull])
		lhs <- rbind(lhs, 1-threshold^vested)
		rhs <- append(rhs, min(1, sum(p.max > threshold)/length(p.max)))
	}
	qr.solve(qr(lhs), rhs)
}

processPvals <- function(pvals)
# Processes p-values to obtain the necessary statistics 
# e.g. maximum p-value, proportion under/over each p-value.
{
    n <- length(pvals)
	stopifnot(n>1L)
    combos <-  t(expand.grid(rep(list(0:1), n)))[,-1,drop=FALSE]
	ng <- length(pvals[[1]])
	stopifnot(all(ng==sapply(pvals, FUN=length)))

    # Compute alpha's and m's.
	alpha <- do.call(pmax, pvals)
	m <- sapply(pvals, FUN=function(x) {
		index <- findInterval(alpha, sort(x))
		pmin(1-alpha, ifelse(index==ng, 1/ng, 1-index/ng))
	})
	return(list(design=combos, p.max=alpha, m=m))
}

