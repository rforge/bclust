\name{teethplot}
\alias{teethplot}
\title{
produces teeth plot useful for demonstating a grouping on clustered subjects}
\description{
The function alone is useless but can be attached to a horizontal dendrogram, a profile plot, or an image plot to show a specified partitioning.
}
\usage{
teethplot(x, teeth.space = 0.25, teeth.lwd = 1)
}
\arguments{
  \item{x}{A \code{bclustvs} object.}
  \item{teeth.space}{The space between two teeth, a value between 0 and 0.25.}
  \item{teeth.lwd}{The thickness of the lines used to draw the teeth.}
}
\details{
The teeth plot for the moment shows the grouping vertically. 
}
\seealso{
\link{profileplot}, \link{ditplot}, \link{dptplot}.
}
\examples{
data(gaelle)
gaelle.id<-rep(1:14,c(3,rep(4,13))) # the first 3 rows are replication of ColWT , the other mutants each has 4 replications.
gaelle.lab<-c("ColWT","d172","d263","isa2","sex4","dpe2","mex1",
"sex3","pgm","sex1","WsWT","tpt","RLDWT","ke103")
gaelle.bclust<-bclust(gaelle,rep.id=gaelle.id,labels=gaelle.lab,
transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68),var.select=TRUE)
#start plotting
layout(matrix(c(1,2),1,2,byrow=TRUE), c(9,1),10, respect=TRUE) # divide plot space into two unequal parts
par(mar=c(0,0,0,2)) # preserve some space for labels in dendrogram plot
plot(as.dendrogram(gaelle.bclust),horiz=TRUE,yaxs="i") #plot the dendrogram
abline(v=gaelle.bclust$cut) # show the optimal allocation with a line on the dendrogram
par(mar=c(0,0,0,0)) # we need no space for teeth plot
teethplot(gaelle.bclust) #show the optimal allocation using teeth plot
}