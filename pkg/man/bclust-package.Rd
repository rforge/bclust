\name{bclust}
\alias{bclust}
\docType{package}
\title{
The package does Bayesian agglomerative clustering with variable selection using spike and slab models in two levels useful for clustering continuous high-dimensional data.
}
\description{
A linear additive model with a diappearing random effect is assumed for data and hierarchical clustering is implemented using the log posterior as a natural distance imposed by the model yielding dendorgrams with probabilistic interpretation. Similicity of model provides estimation of model parameters using empirical Bayes and hence provides a fully automatic model-based approch. The package provides possibility of two distribution families for random effects, Gaussian or asymmetric Laplace. 
}
\details{
\tabular{ll}{
Package: \tab bclust\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-04-20\cr
License: \tab GPL (>=2.0)\cr
LazyLoad: \tab yes\cr
}
bclust, likelihood, bclassify
}
\author{Vahid PARTOVI NIA and Anthony C. DAVISON, Maintainer: Vahid PARTOVI NIA <bclust@probstat.ch>
}
\references{
Fast High-Dimensional Bayesian Classification and Clustering
}
\examples{
data(gaelle)
plot(bclust(gaelle,transformed.par=c(-1.8,-1.7,0.2,0.2,0.6,0))) #transformed.par obtained by maximum likelihood
}
