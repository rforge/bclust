\name{profileplot}
\alias{profileplot}
\title{a plot useful to visualise replicated data}
\description{The profile plot is a flexible function creating profiles of replicated data with many useful options. The resulting plot can be attached to other plots like a horizontal dendrogram plot or a teethplot.}
\usage{
profileplot(x, rep.id, labels = NULL, scale = 1, 
col.names = colnames(x), plot.order = 1:max(rep.id), 
xlab.mar = 5, ylab.mar = 5, xlab.cex = 0.8, 
ylab.cex = 1, profile.col = rep(1, max(rep.id)), 
blob.matrix = matrix(0, ncol = ncol(x), nrow = max(rep.id)), 
blob.col = rev(heat.colors(max(blob.matrix))), blob.cex = 1.5)
}
\arguments{
  \item{x}{
The data matrix with subjects in rows and variables in columns.
}
  \item{rep.id}{
A positive integer vector referring to replication of a subject such that the total number of subjects is \code{max(rep.id)}.
}
  \item{labels}{The labels of each subject. The first element corresponds to the subject that takes value 1 in \code{rep.id}, the second element is the one that takes value 2 in \code{rep.id} 
  and so on. 
}
  \item{scale}{A positive value. If the profile plots are noisy or flat, adjust the scaling factor.}
  \item{col.names}{The data column labels (variable labels) in a string vector.
}
  \item{plot.order}{A positive integer vector, the order that subjects should be plotted. This is useful if you like to attach this plot to a dendrogram or a teeth plot.
}
  \item{xlab.mar}{A positive value, the margin reserved for the column (variable) labels.
}
  \item{ylab.mar}{A positive value, the margin reserved for the subject labels.
}
  \item{xlab.cex}{A positive value, the magnitude of the variable labels.
}
  \item{ylab.cex}{A positive value, the magnitude of the subject labels.
}
  \item{profile.col}{Colours that you may used to plot each profile.
}
  \item{blob.matrix}{An integer matrix denoting where to plot a colour blob and which colour should be used in which location. 
  The number of rows of this matrix should be max( rep.id ) and the number of columns should be \code{ncol(x)}.
  The value 0 plots no blobs, value 1 refers to plotting a blob with a 
  colour specified in \code{blob.col}.
}
  \item{blob.col}{The colours that is used for each value of \code{blob.matrix}}
  \item{blob.cex}{The magnitude of blobs.
}
}
\details{
Some of the options may be useless for just a profile plot, but all of them are beneficial for a good visual representation of data when the plot is attached to a dendrogram or a teethplot.
}
\seealso{
\link{teethplot}, \link{ditplot}, \link{dptplot}.
}
\examples{
data(gaelle)
#take a subset of gaelle data
subgaelle<-gaelle[1:11,]

#use thresholds to define blob colors
blob<-matrix(0,nrow(subgaelle),ncol(subgaelle))
blob[abs(subgaelle)<=0.74]<-0 
blob[abs(subgaelle)>0.74]<-1
blob[abs(subgaelle)>0.94]<-2
blob[abs(subgaelle)>1.28]<-3
profileplot(subgaelle,1:nrow(subgaelle),scale=10,blob.matrix=blob)

#make a profile plot with colored blobs
####### attach a profileplot to a teeth plot
subgaelle.bclust<-bclust(subgaelle,
transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
# divide plot space into two unequal parts
layout(matrix(c(1,2),1,2,byrow=TRUE), c(9,1),10, respect=TRUE) 
# associate a color to each subject
subgaelle.color<-c(rep("blue",3),rep("green",4),rep("magenta",4))
# find appropriate order to plot subjects
subgaelle.order<-order.dendrogram(as.dendrogram(subgaelle.bclust))
#leave some space for labels 
x.mar<-7
y.mar<-5
profileplot(subgaelle,rep.id=1:nrow(subgaelle),scale=10,
profile.col=subgaelle.color,plot.order=subgaelle.order,
xlab.mar=x.mar,ylab.mar=y.mar)
par(mar=c(x.mar,0,0,0))
teethplot(subgaelle.bclust)

####### nowattach it to a dendrogram plot
layout(matrix(c(2,1),1,2,byrow=TRUE), c(2,8),10, respect=TRUE) 
profileplot(subgaelle,rep.id=1:nrow(subgaelle),scale=10,
profile.col=subgaelle.color,plot.order=subgaelle.order,
xlab.mar=x.mar,ylab.mar=y.mar)
par(mar=c(x.mar,0,0,0))
plot(as.dendrogram(subgaelle.bclust),horiz=TRUE,yaxs="i",
axes=FALSE,leaflab="none")
abline(v=subgaelle.bclust$cut)
}
