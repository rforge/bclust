\name{viplot}
\alias{viplot}
\title{variable importance plot}
\description{This function plots varible importances using a barplot.}
\usage{
viplot(varimp, xlab, xlab.mar = 5, ylab.mar = 4, xlab.srt = 90, 
xlab.cex = 1, sort = FALSE, ...)
}
\arguments{
  \item{varimp}{A numeric vector, variable importances.}
  \item{xlab}{A vector of strings. Labels to be plotted on x-axis.}
  \item{xlab.mar}{A positive value. The margin reserved for x-axis.}
  \item{ylab.mar}{A positive value, The margin reserved for y-axis.}
  \item{xlab.srt}{A numeric value, amount of rotation of the x-axis labels.}
  \item{xlab.cex}{A positive value, magnitude of x labels.}
  \item{sort}{Logical, if sorted variables should be plotted.}
  \item{\dots}{The barplot options.}
}
\details{It is not strightforward to rotate x labels of a barplot. This function does it easily.}
\seealso{\link{imp}.}
\examples{
data(gaelle)
gaelle.bclust<-bclust(gaelle,transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
gaelle.imp<-imp(gaelle.bclust)
viplot(varimp=gaelle.imp$var,xlab=imp(gaelle.bclust)$labels)# solid plot

viplot(varimp=gaelle.imp$var,xlab=imp(gaelle.bclust)$labels,
sort=TRUE,col=heat.colors(length(gaelle.imp$var))) # sorted plot
}
