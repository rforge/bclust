\name{meancss}
\alias{meancss}
\title{computes data statistics necessary for evaluation of the model log likelihood}
\description{The function is useful for deriving the maximum likelihood estimates of the model parameters.}
\usage{meancss(x,x.id=1:nrow(x))}
 \arguments{
 \item{x}{The data matrix, clustering individuals in rows, variables in columns.}
 \item{x.id}{A vector of positive integers showing replication of types. The replicates of the same type have the same integers. Apriori considered that data are unreplicated.}
  	}
\value{
\item{mean}{The mean of individuals.}
\item{css}{The corrected sum of squares of individuals.}
\item{repno}{The vector containing the number of replications of individuals according to \code{x.id}.}
}
\details{This function facilitates usage of \code{loglikelihood} function.}
\seealso{
		\link{loglikelihood}
			}
\examples{
data(gaelle)
gaelle.id<-rep(1:14,c(3,rep(4,13))) # first 3 rows replication of ColWT , the other mutants each have four replications
mc.gaelle<-meancss(gaelle,gaelle.id)
loglikelihood(x.mean=mc.gaelle$mean,x.css=mc.gaelle$css,
repno=mc.gaelle$repno,transformed.par=rep(0,6))
}

