\name{giut}
\alias{giut}

\title{Generalized IUT}
\description{Perform a generalized intersection-union test for multiple sets of p-values}

\usage{
giut(..., threshold=0.5)
}

\arguments{
\item{...}{two or more numeric vectors of p-values of the same length}
\item{threshold}{a numeric scalar indicating the p-value threshold with which to estimate proportions}
}

\details{
The generalized IUT is a heuristic method that computes (probably) conservative
p-values against the union null hypothesis. Consider each configuration of
true/false nulls across comparisons, and the proportion of genes corresponding
to each configuration. The generalized approach will search across the
proportion space to identify the locally maximum p-value.

The landscape is fairly bumpy so convergence is not guaranteed. Rather,
optimization is initiated from the most relevant point, i.e., estimates of the
proportions obtained using a multi-comparison extension of Storey's method with
lambda set at \code{threshold}. Maximization will then, hopefully, give a
p-value above the true value.

This method depends on a large number of tests to obtain precise estimates
of the starting proportions, as well as the values for \code{m1}, \code{m2}, 
and so on for each comparison (as described for \code{\link{diut}}). It also
assumes that the alternative distribution is the same within each comparison.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{diut}}
}

\examples{
p1 <- runif(10000)
p2 <- runif(10000)
p3 <- runif(10000)
head(giut(p1, p2, p3))
     
p1[1:1000] <- 1e-8
head(giut(p1, p2, p3))
   
p2[5001:10000] <- rbeta(5000, 1, 20)
head(giut(p1, p2, p3))
}
