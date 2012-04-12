\name{dptplot}
\alias{dptplot}
\title{dendrogram-profile-teeth plot}
\description{This is a plot function suitable for visualisation of a \code{bclustvs} object when the \code{bclust} function is used on replicated data, but can be applied for unreplicated data too. 
The function attaches a coloured horizontal dendrogram to the left side of a profile plot with the optimal grouping highlighted by a teethplot on the right.}
\usage{
dptplot(x, scale = 1, xlab = colnames(x$data), ylab = x$labels, 
xlab.cex = 1, ylab.cex = 1, dendrogram.lwd = 1, 
dendrogram.size = 2, xlab.mar = 3, ylab.mar = 3, 
horizbar.plot = FALSE, 
horizbar.col = rev(c(heat.colors(5)[-4], "white")), 
horizbar.distance = 4,  varimp = rep(0, ncol(x$data)), 
horizbar.size = 0.5, vertbar = NULL, 
vertbar.col = rainbow(max(vertbar)), 
teeth.size = 0.25, plot.width = 10)}
\arguments{
  \item{x}{A \code{bclustvs} object.}
  \item{scale}{A positive value. If the profile plots are too noisy or too flat, adjust the scaling factor.}
  \item{xlab}{A string vector, the labels for the clustering types.}
  \item{ylab}{A string vector, the variable labels.}
  \item{xlab.cex}{A positive value, the magnitude of the type labels.}
  \item{ylab.cex}{A positive value, the magnitude of the variable labels.}
  \item{dendrogram.lwd}{A positive value, the thickness of lines used to plot the dendrogram.}
  \item{dendrogram.size}{A positive value, the size of the dendrogram plot.}
  \item{xlab.mar}{A positive value, the margin reserved to write variable labels.}
  \item{ylab.mar}{A positive value, the margin reserved to write type labels.}
\item{horizbar.plot}{A logical value. If \code{TRUE} a horizontal bar is plotted according to the \code{varimp}.}
  \item{horizbar.col}{Colours of the horizontal bar.}
  \item{horizbar.distance}{A positive value, the distance between the horizontal bar and the profile plot.}
    \item{varimp}{A numerical vector denoting the importance of variables. We propose to use \code{imp} function to get the values. If it is specified, the variables will be ordered respect to this vector.}
  \item{horizbar.size}{A positive value, the size of the horizontal bar.}
  \item{vertbar}{A positive integer vector that may be used to draw an additional vertical bar on the right of the teethplot. This may be helpful to represent another optional grouping on the data.}
  \item{vertbar.col}{The colours used to plot the additional vertical bar.}
  \item{teeth.size}{A positive integer, the size of the teeth plot.}
  \item{plot.width}{A positive integer, the width of the whole plot. If the plot region is unbalanced in width and height, adjust this value.}
}
\details{The \code{varimp} is assumed to be the log Bayes factor and therefore categorised according to Kass and Raftery (1995) for a better visualisation. If \code{varimp} is not 
result of the \code{imp} function, keep \code{horizbar.plot=FALSE}.}
\references{
Kass and Raftery (1995) Bayes Factors, Journal of the American Statistical Association, Vol. 90, pp. 773--795.
}
\seealso{
\link{ditplot}, \link{profileplot}, \link{teethplot}, \link{viplot}.
}
\examples{
data(gaelle)

gaelle.id<-rep(1:14,c(3,rep(4,13))) 
# first 3 rows replication of ColWT, 4 for the rest

gaelle.lab<-c("ColWT","d172","d263","isa2",
"sex4","dpe2","mex1","sex3","pgm","sex1","WsWT","tpt","RLDWT","ke103")
gaelle.bclust<-bclust(gaelle,rep.id=gaelle.id,labels=gaelle.lab,
transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
dptplot(gaelle.bclust,scale=5,varimp=imp(gaelle.bclust)$var,
horizbar.plot=TRUE)
#replicated clustering


gaelle.bclust<-bclust(gaelle,
transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
dptplot(gaelle.bclust,scale=10,varimp=imp(gaelle.bclust)$var,
horizbar.plot=TRUE,plot.width=5,horizbar.size=0.2,ylab.mar=4)
#unreplicated clustering

wildtype<-rep(1,55) #initiate a vector
wildtype[c(1:3,48:51,40:43)]<-2 #associate 2 to wildtypes
dptplot(gaelle.bclust,scale=10,varimp=imp(gaelle.bclust)$var,
horizbar.plot=TRUE,plot.width=5,horizbar.size=0.2,vertbar=wildtype,
vertbar.col=c("white","violet"),ylab.mar=4)

}
