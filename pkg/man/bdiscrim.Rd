\name{bdiscrim}
\alias{bdiscrim}
\title{discrimination using a Bayesian linear model}
\description{This function provides a discrimination method equivalent to the \code{bclust} clustering function.}
\usage{bdiscrim(training, training.id, 
training.labels = NULL, predict, 
predict.label = rownames(predict)[1], 
effect.family = "gaussian", var.select = TRUE, 
transformed.par, priorprob = rep(1, max(training.id) + 1))}
\arguments{
  \item{training}{A numeric matrix of the training set with observations in rows and variables in columns.}
  \item{training.id}{A vector consisting of positive integer elements having the
 same length as the number of rows of \code{training}. This vector identifies observations of a class such that the total number of classes is \code{max(training.id)}. }
  \item{training.labels}{A vector of strings
 referring to the labels of each class. The length of the
 vector should match to \code{max(training.id)}. The first element corresponds to the label of the type having the smallest integer value in \code{training.id}, the second element refers to 
 the label of the type having the second smallest integer in \code{training.id}, and so on.} 
  \item{predict}{A single discriminating type to be classified to one of the classes, with replications of the discriminating type in rows and variables in columns. The number of variables should be the same as in \code{training}. If it is not specified just the variable and variable-class importances will be calculated.}
  \item{predict.label}{A single string, the label of the discriminating type.}
  \item{effect.family}{Distribution family of the disappearing random components of the model. Choices are "gaussian" or "alaplace"
 allowing Gaussian or asymmetric Laplace family, respectively.}
  \item{var.select}{A logical value, \code{TRUE} for fitting models that define
   a spike-and-slab distribution in variable level and allows Bayesian variable selection.}
  \item{transformed.par}{The transformed
 model parameters in a vector. The length of the vector depends on the
 chosen model and the availability of variable selection. The log transformation is supposed to be applied for the variance parameters, the identity for the mean, and the logit for the proportions. The function \link{loglikelihood} can be used to estimate them from the data.}
  \item{priorprob}{The prior probabilities for each class in a vector. The length of the vector should be the number of classes \code{max(training.id)} plus one (the prior probability that the predict set
  raises its own class).  If nothing is specified a uniform discrete prior is considered.}
}
\details{
The function calls internal \code{C} functions depending on the chosen
model. The C-stack of the system may overflow if you have a large
dataset. You may need to adjust the stack before running \code{R} using your
operational system command line. If you use Linux open console type
\code{>ulimit -s unlimited}, then run \code{R} in the same console. 

We assumed a Bayesian linear model for classification being \deqn{y_{vctr}=\mu+\eta_{vct}+\delta_v \gamma_{vc}\theta_{vc}+\varepsilon_{vctr}}{y_{vctr}=m+h_{vct}+d_{v}*g_{vc}*t_{vc}+e_{vctr}} where \eqn{y_{vctr}} is the available data on variable \eqn{v}, 
cluster \eqn{c}, type \eqn{t}, and replicate \eqn{r}; \eqn{\eta_{vct}}{h_{vct}}
is the between-type error, \eqn{\theta_{vc}}{t_{vc}} is the disappearing random component controlled by the Bernoulli variables \eqn{\delta_v}{d_{v}} with success probability \eqn{q} and \eqn{\gamma_{vc}}{g_{vc}} with 
success probability \eqn{p}; and \eqn{\varepsilon_{vctr}}{e_{vctr}} is the between-replicate error. The function computes the posterior probability that the \code{predict} data share the same \eqn{\theta_{vc}}. 
This function also consideres that the \code{predict} data may arise its own class.}
\value{
\item{probs}{The posterior probabilities in a matrix with one row.}
\item{var}{The variable importances, each being a log Bayes factor.}
\item{varclass}{The variable-class importances, each being a log Bayes factor.}
}
\seealso{\link{loglikelihood}, \link{bclust}, \link{profileplot}.}
\examples{
data(gaelle)
gaelle.id<-rep(1:13,rep(4,13)) # all mutants have 4 replicates
gaelle.lab<-c("d172","d263","isa2",
"sex4","dpe2","mex1","sex3","pgm","sex1","WsWT","tpt","RLDWT","ke103")

gaelle.bdiscrim<-bdiscrim(gaelle[-(1:3),],
training.id=gaelle.id,training.labels=gaelle.lab,
transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68),
predict=gaelle[1:3,], predict.label="ColWT")
# classify ColWT to one of the types

par(mfrow=c(1,1)) #retreive graphic defaults 
par(mar = c(5, 4, 4, 2) + 0.1) # leave some space for labels

ColWT.prob<-as.vector(gaelle.bdiscrim$probs)

# plots discrimination probabilities
bp <- barplot(ColWT.prob,ylim=c(0,1)) #plot bars
title("ColWT Discrimination")
text(bp, par("usr")[3]-0.05, srt = 90,adj=1,
labels = colnames(gaelle.bdiscrim$probs), xpd = TRUE) 
#plot variable labels

mtext(1, text = "Mutant", line = 4,cex=1.5)# add x axis label
mtext(2, text = "Probability", line = 3,cex=1.2)# add y axis labels
abline(h=1/length(ColWT.prob)) 
# draw plot a as a reference line prior probabilities line



# plots sorted discrimination probabilities
par(mfrow=c(1,1)) #retreive graphic defaults 
par(mar = c(5, 4, 4, 2) + 0.1) # leave some space for labels
bp <- barplot(sort(ColWT.prob,decreasing=TRUE),ylim=c(0,1),
col=heat.colors(length(ColWT.prob))) #plot bars
text(bp, par("usr")[3]-0.05, srt = 90,adj=1,
labels = colnames(gaelle.bdiscrim$probs)[order(ColWT.prob,decreasing=TRUE)], 
xpd = TRUE) #plot variable labels
mtext(1, text = "Mutant", line = 4,cex=1.5)# add x axis label
mtext(2, text = "Probability", line = 3,cex=1.2)# add y axis labels
abline(h=1/length(ColWT.prob))


varclassimp<-gaelle.bdiscrim$varclass
#use thresholds to define blob colors

     blob<-matrix(0,nrow(varclassimp),ncol(varclassimp))
     blob[varclassimp<=0]<-0 
     blob[varclassimp>0]<-1
     blob[varclassimp>1]<-2
     blob[varclassimp>3]<-3
     blob[varclassimp>5]<-4
#log bayes factor thresholding 

varimp<-gaelle.bdiscrim$var
varcol<-rep(0,ncol(gaelle))
varcol[varimp>0]<-1


var.order<-order(gaelle.bdiscrim$var,decreasing=TRUE)
profileplot(x=gaelle[-(1:3),var.order],rep.id=gaelle.id,
labels=gaelle.lab,scale=10,blob.matrix=blob,ylab.mar=5,xlab.mar=7)
#plot var class importance on profile plot using blobs
# and sort variables according to variable importance values



viplot(varimp=varimp, xlab=colnames(gaelle)[var.order],
xlab.mar=10,sort=TRUE,col=varcol[var.order])
#plot sorted variable importances 
}
