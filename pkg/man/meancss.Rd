\name{meancss}
\alias{meancss}
\title{computes statistics necessary for the evaluation of the log likelihood}
\description{The function is useful for deriving arguments of the \link{loglikelihood} function.}
\usage{meancss(x,rep.id=1:nrow(x))}
 \arguments{
 \item{x}{The data matrix, subjects in rows, variables in columns.}
 \item{rep.id}{A vector of positive integers refering to replication of the types. The same integer is associated to the the replicates of the same type. 
 Not specifying this vector preproposes that the data are unreplicated.}
  	}
\value{
\item{mean}{The mean of types.}
\item{css}{The corrected sum of squares of types.}
\item{repno}{The vector containing the number of replications of types according to \code{rep.id}.}
}
\details{This function facilitates the usage of the \code{loglikelihood} function.}
\seealso{
		\link{loglikelihood}.
			}
\examples{
data(gaelle)
gaelle.id<-rep(1:14,c(3,rep(4,13))) # first 3 rows replication of ColWT , the other mutants each have four replications
mc.gaelle<-meancss(gaelle,gaelle.id)
loglikelihood(x.mean=mc.gaelle$mean,x.css=mc.gaelle$css,
repno=mc.gaelle$repno,transformed.par=rep(0,6)) # evalutes likelihood at rep(0,6)
}

