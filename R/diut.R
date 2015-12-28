diut <- function(p1, p2) 
# A special case, when you have two sets of p-values. The 'd' is for 'di-',
# obviously.  Protection is provided so that it false back on the maximum
# p-value if it gets too large.
{
	all.p <- list(p1, p2)
	out <- processPvals(all.p)
	out$p.max * pmin(1, 2 - rowSums(out$m) - out$p.max)
}

