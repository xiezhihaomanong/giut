require(giut)
computePval <- function(alpha, prop, design, M)
{
	nullprop <- as.vector(design %*% prop)
	betas <- (M - (1-alpha)*nullprop)/(1-nullprop)
	mult <- (1-betas)*(1-design) + alpha*design
	log(sum( exp(colSums(log(mult))) * prop )/sum(prop))
}

tester <- 62L
ref.iut <- function(..., threshold=0.5)
# This function is otherwise identical to the other one, except 
# that you can change 'tester' to see the various workings out
# for a particular analysis.
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
	
		if (o[i]==tester){ 
			print(sprintf("Gene: %i (%i) %.10f", o[i], i, cur.alpha));
		}	
			
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
			if (o[i]==tester) {
				print(paste(sprintf("%.10f", theta), collapse=" ")) 
			}
			theta <- constrOptim(theta, ui=ui, ci=cur.ci, f=function(prop) {
				diff0 <- theta0 - prop
				sum(diff0 * diff0)
			}, method="Nelder", control=list(reltol=1e-8,trace=o[i]==tester))$par
			past <- theta
		}
		if (o[i]==tester) {
			print(paste(sprintf("%.10f", theta), collapse=" ")) 
		}

		# Maximizing it.
		soln <- constrOptim(theta, ui=ui, ci=cur.ci, f=computePval, method="Nelder",
			alpha=cur.alpha, design=out$design, M=cur.m, control=list(fnscale=-1, reltol=1e-8))#trace=o[i]==tester))
		pval[o[i]] <- exp(soln$value)
		if (o[i]==tester) { 
			print(pval[o[i]]) 
			print(paste(sprintf("%.10f", soln$par), collapse=" ")) 
		}
	}
	return(pval)
}

set.seed(1231)
ngenes <- 100
p1 <- rbeta(ngenes, 1, 20)
p2 <- runif(ngenes)
p2[81:100] <- 0
p3 <- runif(ngenes)
p3[51:100] <- rbeta(10, 1, 20)

ref <- ref.iut(p1, p2, p3)
test <- giut(p1, p2, p3)

order(abs(ref-test), decreasing=TRUE)
order(abs(pmax(p1, p2, p3)), decreasing=TRUE)

# For the chosen analysis, the R and C++ implementations diverge during the
# least-squares minimization! Of all places. Starting inputs are the same,
# and starting barrier function values are the same. Divergence in the middle 
# of the trace that it's a difference in the underlying compiled implementation,
# somehow...
