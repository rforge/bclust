\name{bclustvs}
\alias{bclustvs}
\title{bclustvs (Bayesian CLUSTering with Variable Selection) is a class}
\description{The \code{bclustvs} object can be regarded as a Bayesian extension of the \code{hclust} object that additionally contains information of a Bayesian dendrogram. 
You can convert any \code{bclustvs} object to a \code{dendrogram} object using \code{as.dendrogram}.}
\value{
	\item{data}{The data matrix, reordered according to \code{rep.id}.}
	\item{repno}{The number of replicates of the values of \code{rep.id}}
	\item{merge}{The merge matrix in \code{hclust} object format.}
	\item{height}{A monotone vector referring to the height of the
	constructed tree.}
	\item{logposterior}{The log posterior for each merge.}
	\item{clust.number}{The number of clusters in each merge.}
	\item{cut}{The value of the height corresponding to the maximum of the log posterior in agglomerative path.}
	\item{transformed.par}{The transformed values of the model parameters. 
	  The log transformation is applied for the variance parameters, the identity for the mean, and the logit for the proportions.}
	\item{labels}{The labels associated to each clustering type.}
	\item{effect.family}{The distribution assigned to the disappearing random effect.}
	\item{var.select}{Availability of variable selection chosen in the function arguments.}
      }

