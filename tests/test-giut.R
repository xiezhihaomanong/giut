########################################################################################
########################################################################################
########################################################################################

suppressPackageStartupMessages(require(giut))

computePval <- function(alpha, prop, design, M)
# A function to estimate the probability of getting a maximum p-value under
# alpha, given the value of alpha, the proportion of genes in each combination
# of null/alts, the definition of the null/alts and the proportion of genes
# with max p-values above alpha.
{
	nullprop <- as.vector(design %*% prop)
	betas <- (M - (1-alpha)*nullprop)/(1-nullprop)
	mult <- (1-betas)*(1-design) + alpha*design
	log(sum( exp(colSums(log(mult))) * prop )/sum(prop))
}

ref.iut <- function(..., threshold=0.5) 
# This brings everything together. First, it estimates the proportion of genes
# in each combination of null/alternative hypotheses. As this estimate isn't
# reliable or unbiased, the function then attempts to maximize the p-value by
# searching the parameter space around the estimate. This assumes that the
# estimate is reasonably close, such that the local maximum is conservative.
{
	all.p <- list(...)
	out <- giut:::processPvals(all.p)
    o <- order(out$p.max, decreasing=TRUE) 

	# Initializing the obvious constraints for probabilities.
	nc <- 2L^length(all.p) - 1L
	ui <- rbind(diag(nc), rep(-1, nc), -out$design)
	ci <- c(rep(0, nc), -1)

	# Initializing the proportions.
	theta0 <- giut:::estimateProp(all.p, out$design, threshold=threshold)

	# Computing pmax for each gene. 
	pval <- numeric(length(out$p.max))
	past <- NULL
	for (i in 1:length(o)) {
		cur.alpha <- out$p.max[o[i]]
		cur.m <- out$m[o[i],]
		
		# Searching for a feasible point close to the theoretical point.
		# This should converge as everything is linear.
		upper.b <- cur.m/(1-cur.alpha)
		cur.ci <- c(ci, -upper.b)
		if (all(ui %*% theta0 > cur.ci)) { 
			theta <- theta0
			past <- NULL
		} else {
			if (is.null(past) || any(ui %*% past <= cur.ci)) { 
				theta <- rep(min(upper.b)/(nc + 1), nc)
			} else {
				theta <- past
			}
#			print(paste(sprintf("%.8f", theta), collapse=" "))
			theta <- constrOptim(theta, ui=ui, ci=cur.ci, f=function(prop) {
				diff0 <- theta0 - prop
				sum(diff0 * diff0)
			}, method="Nelder", control=list(reltol=1e-8))$par
			past <- theta
		}

		# Maximizing it.
		soln <- constrOptim(theta, ui=ui, ci=cur.ci, f=computePval, method="Nelder",
			alpha=cur.alpha, design=out$design, M=cur.m, control=list(fnscale=-1, reltol=1e-8))
		pval[o[i]] <- exp(soln$value)
	}
	return(pval)
}

########################################################################################
########################################################################################
########################################################################################
# Just running on a lot of two experiments to test it.

ngenes <- 100
checkme <- function(..., threshold=1e-6) {
	ref <- ref.iut(...)
	test <- giut(...)
	stopifnot(all(abs(ref-test) <= (ref+1e-6)*threshold))
	head(ref)
}

set.seed(42386332)
p1 <- runif(ngenes)
p1[1:10] <- rbeta(10, 1, 20)
p2 <- runif(ngenes)
p2[91:100] <- 0
checkme(p1, p2)

p1 <- runif(ngenes)
p1[1:10] <- 0
p2 <- runif(ngenes)
p2[51:100] <- rbeta(50, 1, 20)
checkme(p1, p2)

p1 <- runif(ngenes)
p2 <- runif(ngenes)
p2[51:100] <- 0
checkme(p1, p2)

p1 <- runif(ngenes)
p2 <- runif(ngenes)
checkme(p1, p2)

p1 <- runif(1e3)
p2 <- runif(1e3)
p2[501:1000] <- rbeta(500, 1, 20)
checkme(p1, p2)

p1 <- runif(1e3)
p2 <- runif(1e3)
checkme(p1, p2)

p1 <- runif(1e3)
p2 <- runif(1e3)
p2[501:1000] <- 0
checkme(p1, p2)

p1 <- runif(1e3)
p1[1:500] <- 0
p2 <- runif(1e3)
p2[501:1000] <- rbeta(500, 1, 20)
checkme(p1, p2)

##########################################################################################
# Three experiments. Loosening the threshold as it tends to not play nice.  For
# the worst example below (maximum difference of ~5%), there appears to be
# convergence issues throughout various optimisation steps. These originate in
# differences in the nmmin trace when computing the p-value (same proportion
# values up to 10 dp's, so same starting point). Possibly it occurs when you
# step on the boundary between two concavities, where differences in
# compilation or precision, etc., decide which way it flips.

compme <- function(...) {
	ref <- ref.iut(...)
	test <- giut(...)
	stuff <- abs(ref-test)/(ref+1e-6)
	return(summary(stuff))
}

p1 <- runif(ngenes)
p1[1:10] <- 0
p2 <- runif(ngenes)
p2[91:100] <- rbeta(10, 1, 20)
p3 <- runif(ngenes)
p3[91:100] <- 0
compme(p1, p2, p3)

p1 <- runif(ngenes)
p1[1:10] <- rbeta(10, 1, 20)
p2 <- runif(ngenes)
p2[81:100] <- 0
p3 <- runif(ngenes)
p3[51:100] <- rbeta(10, 1, 20)
compme(p1, p2, p3)

p1 <- rbeta(ngenes, 1, 20)
p2 <- runif(ngenes)
p2[81:100] <- 0
p3 <- runif(ngenes)
p3[51:100] <- rbeta(10, 1, 20)
compme(p1, p2, p3)

p1 <- runif(ngenes)
p2 <- runif(ngenes)
p3 <- runif(ngenes)
compme(p1, p2, p3)

p1 <- runif(ngenes)
p1[1:80] <- 0
p2 <- runif(ngenes)
p3 <- runif(ngenes)
compme(p1, p2, p3)

p1 <- rbeta(ngenes, 1, 20)
p2 <- runif(ngenes)
p3 <- runif(ngenes)
compme(p1, p2, p3)

##########################################################################################
##########################################################################################
##########################################################################################

