\name{imp}
\alias{imp}
\title{
calculates variable and variable-cluster importances
}
\description{
The function computes the log Bayes factors for the hypothesis H0: the variable or the variable-cluster combination is useful for clustering  against H1: the variable or the variable-cluster combination is useless. The Bayes factors are computed for the optimal allocation found by the \code{bclust} function.
}
\usage{
imp(x)
}
\arguments{
  \item{x}{ A \code{bclustvs} object.}
}
\value{
\item{var}{A vector being the log Bayes factor of \eqn{\delta_v=1}{d_{v}=1} against \eqn{\delta_v=0}{d_{v}=0}, see \link{bclust} for details.}
\item{varclust}{A vector being the log Bayes factor of \eqn{\gamma_{vc}=1}{g_{vc}=1} against \eqn{\gamma_{vc}=0}{g_{vc}=0}, see \link{bclust} for details.}
\item{repno}{The number of replicates producing each row of \code{varclust}.}
\item{labels}{The vector of variable labels extracted from the \code{bclustvs} object.}
\item{order}{The order of \code{var} useful to sort \code{var}, \code{varclust}, and \code{labels}.}
}
\seealso{\link{bclust}.}
\examples{
data(gaelle)

gaelle.id<-rep(1:14,c(3,rep(4,13))) # first 3 rows replication of ColWT , the other mutants each 

gaelle.bclust<-bclust(gaelle,rep.id=gaelle.id,transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68),var.select=TRUE)

gaelle.imp<-imp(gaelle.bclust)

#plot the variable importances 
par(mfrow=c(1,1)) #retreive graphic defaults 

mycolor<-gaelle.imp$var
mycolor<-c()
mycolor[gaelle.imp$var>0]<-"black"
mycolor[gaelle.imp$var<=0]<-"white"

viplot(var=gaelle.imp$var,xlab=gaelle.imp$labels,col=mycolor)
#plot important variables with balck

viplot(var=gaelle.imp$var,xlab=gaelle.imp$labels,
sort=TRUE,col=heat.colors(length(gaelle.imp$var)),
xlab.mar=10,ylab.mar=4)
mtext(1, text = "Metabolite", line = 7,cex=1.5)# add x axis label
mtext(2, text = "Log Bayes Factor", line = 3,cex=1.2)# add y axis labels
#sort importnaces and use heat colors add some labels to the x and y axes
}
