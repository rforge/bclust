\name{imp}
\Rdversion{1.1}
\alias{imp}
\title{
calculates variable and variable-cluster importances
}
\description{
The function computes the log Bayes factors for the hypothesis H0: variable or variable-cluster combination is useful for clustering  according to 
for the optimal allocation found by agglomerative method.
}
\usage{
imp(x)
}
\arguments{
  \item{x}{ A \code{bclustvs} object created by \code{bclust} function.}
}
\value{
\item{var}{A sorted vector being the log Bayes factor of \eqn{\delta}=1 against \eqn{\delta}=0 , see \link{bclust} for details.}
\item{varclust}{A vector being the log Bayes factor of \eqn{\gamma}=1 against \eqn{\gamma}=0, ordered according to \code{var}, see \link{bclust} for details.}
\item{labels}{A vector of variable labels ordered according to \code{var}.}
\item{order}{The order of \code{var} used to sort \code{varclust} and \code{labels}.}
}
\seealso{bclust}
\examples{
data(gaelle)

gaelle.id<-rep(1:14,c(3,rep(4,13))) # first 3 rows replication of ColWT , the other mutants each 

gaelle.bclust<-bclust(gaelle,rep.id=gaelle.id,transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68),var.select=TRUE)

gaelle.imp<-imp(gaelle.bclust)


par(mar = c(10, 4, 4, 2) + 0.1) # leave some space for labels
bp <- barplot(gaelle.imp$var) #plot bars
text(bp, par("usr")[3] - 0.5, srt = 90, adj = 1,
labels = gaelle.imp$labels, xpd = TRUE) #plot variable labels
mtext(1, text = "Metabolite", line = 7,cex=1.5)# add x axis label
mtext(2, text = "Log Bayes Factor", line = 3,cex=1.2)# add y axis labels

#now bars are ordered respect to variable importance
par(mar = c(10, 4, 4, 2) + 0.1) # leave some space for labels
bp <- barplot(gaelle.imp$var[gaelle.imp$order],col=heat.colors(length(gaelle.imp$var))) #plot color bars
text(bp, par("usr")[3] - 0.5, srt = 90, adj = 1,
labels = gaelle.imp$labels[gaelle.imp$order], xpd = TRUE)#plot variable labels
mtext(1, text = "Metabolite", line = 7,cex=1.5)# add x axis label
mtext(2, text = "Log Bayes Factor", line = 3,cex=1.2)# add y axis labels
}
