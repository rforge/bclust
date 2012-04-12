\name{ditplot}
\alias{ditplot}
\title{dendrogram-image-teeth plot}
\description{This is a handy function to plot a \code{bclustvs} object.  The function attaches a coloured horizontal dendrogram to the left side of an image plot with the optimal grouping 
highlighted by a teethplot on the right.}
\usage{
ditplot(x, xlab = colnames(x$data), 
ylab = x$labels, xlab.cex = 1, 
ylab.cex = 1, dendrogram.lwd = 1, dendrogram.size = 2, 
xlab.mar = 3, ylab.mar = 3, image.col = rainbow(20), 
horizbar.plot = FALSE, 
horizbar.col = rev(c(heat.colors(5)[-4], "white")), 
horizbar.distance = 4,varimp = rep(0, ncol(x$data)), 
horizbar.size = 0.5, vertbar = NULL, 
vertbar.col = rainbow(max(vertbar)), 
teeth.size = 0.25, plot.width = 10)
}
\arguments{
  \item{x}{A \code{bclustvs} object.}
  \item{xlab}{A vector of strings elements. The labels for the clustering types automatically extracted from \code{x}.}
  \item{ylab}{A vector of strings. The variable labels automatically extracted from \code{x}.}
  \item{xlab.cex}{A positive value, the magnitude of the clustering type labels.}
  \item{ylab.cex}{A positive value, the magnitude of the variable labels.}
  \item{dendrogram.lwd}{A positive value, the thickness of lines used to plot the dendrogram.}
  \item{dendrogram.size}{A positive value, the size of the dendrogram plot.}
  \item{xlab.mar}{A positive value, the margin reserved to write variable labels.}
  \item{ylab.mar}{A positive value, the margin reserved to write type labels.}
  \item{image.col}{Colours used for the image plot.}
  \item{horizbar.plot}{A logical value, if \code{TRUE} a horizontal bar is plotted according to categorised \code{varimp}.}
  \item{horizbar.col}{Colours used for the horizontal bar.}
  \item{horizbar.distance}{A positive value, the distance between the horizbar and the image plot.}
  \item{varimp}{A numerical vector denoting the importance of variables. We propose to use the \code{imp} function to compute these values. If it is specified, the variables will be ordered respect to this vector.}
  \item{horizbar.size}{A positive value, the size of the horizontal bar.}
  \item{vertbar}{A positive integer vector that may be used to draw an additional vertical bar on the right of the teeth plot. This may be helpful to represent another optional grouping on the data.}
  \item{vertbar.col}{The colours used to plot the additional vertical bar.}
  \item{teeth.size}{A positive value, the size of the teeth plot.}
  \item{plot.width}{A positive value, the width of the whole plot. If the plot region is unbalanced in width and height, adjust this value.}
  }
\details{The \code{varimp} is assumed to be log Bayes factors and therefore categorised according to Kass and Raftery (1995) for a better visualisation. }
\references{
Kass and Raftery (1995) Bayes Factors, Journal of the American Statistical Association, Vol. 90, pp. 773--795.
}
\seealso{
\link{dptplot}, \link{teethplot},\link{profileplot},\link{viplot}.
}
\examples{
data(gaelle)

gaelle.bclust<-bclust(gaelle,
transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
ditplot(gaelle.bclust,varimp=imp(gaelle.bclust)$var,horizbar.plot=TRUE,
plot.width=5,horizbar.size=0.2,ylab.mar=4)
#unreplicated clustering

wildtype<-rep(1,55) #initiate a vector
wildtype[c(1:3,48:51,40:43)]<-2 #associate 2 to wildtypes
ditplot(gaelle.bclust,varimp=imp(gaelle.bclust)$var,horizbar.plot=TRUE,
plot.width=5,horizbar.size=0.2,vertbar=wildtype,
vertbar.col=c("white","violet"),ylab.mar=4)
#mark wildtype plants using violet 

gaelle.id<-rep(1:14,c(3,rep(4,13))) 
# first 3 rows replication of ColWT, 4 for the rest
gaelle.lab<-c("ColWT","d172","d263","isa2",
"sex4","dpe2","mex1","sex3","pgm","sex1","WsWT","tpt","RLDWT","ke103")
gaelle.bclust<-bclust(gaelle,rep.id=gaelle.id,
labels=gaelle.lab,transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
ditplot(gaelle.bclust,varimp=imp(gaelle.bclust)$var,horizbar.plot=TRUE)
#replicated clustering
}
