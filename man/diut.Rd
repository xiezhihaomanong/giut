\name{diut}
\alias{diut}

\title{Restricted IUT}
\description{Perform a restricted intersection-union test with two sets of independent p-values}

\usage{
diut(p1, p2)
}

\arguments{
\item{p1, p2}{a numeric vector of p-values of the same length}
}

\details{
This function performs the restricted intersection-union test on two sets of
p-values. The sampling of p-values must be independent between vectors  though
the  configuration of true/false null hypotheses need not be. Mixtures of
alternative p-value distributions are tolerated, as are correlations between
tests within each set.

To calculate the restricted IUT p-value, the maximum of \code{p1} and \code{p2}
is first obtained and defined as \code{pm}. The proportion of p-values that are
greater than this maximum value is determined in each set and denoted as
\code{m1} and \code{m2}. The IUT p-value is then defined as \code{min(pm,
pm*(2-m1-m2-pm))}.

This function will work best with large vectors in order to obtain precise
estimates for \code{m1} and \code{m2}. Nonetheless, the estimated p-value will
be an unbiased estimate of the upper bound of the error rate with respect to
imprecision in \code{m1} and \code{m2}.
}

\author{
Aaron Lun
}

\examples{
p1 <- runif(10000)
p2 <- runif(10000)
head(diut(p1, p2))

p1[1:1000] <- 1e-8
head(diut(p1, p2))

p2[5001:10000] <- rbeta(5000, 1, 20)
head(diut(p1, p2))
}

