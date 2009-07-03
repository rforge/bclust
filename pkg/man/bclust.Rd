\name{bclust}

\alias{bclust}

\title{Bayesian agglomerative clustering for high dimensional data with variable selection.}

\description{The function clusters data saved in a matrix using an additive linear model with disappearing random effects. 
The model has built-in spike-and-slab components which quantifies important variables for clustering.}

\usage{bclust(x,transformed.par,rep.id=1:nrow(x),effect.family="gaussian",var.select=TRUE,labels=NULL)
}


\arguments{
  \item{x}{A numeric matrix of clustering indiviuals in rows
    and variables in columns.}
  \item{transformed.par}{The transformed
 model parameters in a vector. The length of the vector depends on the
 chosen model and the availability of variable selection.}
 
\item{rep.id}{A vector consisting of positive integer elements having the
 same length as the number of rows of \code{x}. This vector identifies replicates of a type such that the total number of clustering types is \code{max(rep.id)}. If nothing is declared the function
 presupposes that data are unreplicated, that is each row of \code{x} is a clustering type.}
\item{effect.family}{Distribution family of the disappearing random components. Choices are "gaussian" or "alaplace"
 allowing Gaussian or asymmetric Laplace family, respectively. }

 \item{var.select}{A logcial value, TRUE for fitting models that define
   spike-and-slab distribution in variable level and allows Bayesian variable selection.}

 \item{labels}{A vector of strings
 refering to the labels of clustering individuals. The length of the
 vector should match to \code{max(rep.id)}. The first element corresponds to the label of the type having the smallest integer value in \code{rep.id}, the second element referes to 
 the label of the type having the second smallest integer in \code{rep.id}, and so on.}} 

\details{The function calls internal C functions depending on the chosen
model. The C-stack of the system may overflow if you have a large
dataset. You may need to adjust the stack before running R using your
operational system command line. If you use Linux open console and type
\code{>ulimit -s unlimited}, then run \code{R} in the same console. 

We assumed Bayesian model for clustering being \deqn{y_{vctr}=\mu+\eta_{vct}+\delta_v \gamma_{vc}\theta_{vc}+\varepsilon_{vctr};} where \eqn{y} is the avilable data on variable \eqn{v}, 
cluster \eqn{c}, clustering type \eqn{t}, and replicate \eqn{r}; \eqn{\eta} 
is between type error, \eqn{\theta} is the disppearing random component controlled by bernoulli variables \eqn{\delta} with succes probability \eqn{q} and \eqn{\gamma} with 
success probability \eqn{p}; and \eqn{\varepsilon} is the between replicate error. The types inside a cluster share the same \eqn{\theta_{vc}}, but may araise with a different  \eqn{\eta_{cvt}}, 
for more details see the package website \url{http://bclust.probstat.ch}.
}

\value{
	\item{data}{The data matrix reordered according to \code{rep.id}.}
	\item{repno}{The number of replicates of the values of \code{rep.id}}
	\item{merge}{The merge matrix in \code{hclust} format.}
	\item{height}{A monotone vector reffering to the height of the
	constructed tree.}
	\item{logposterior}{The log posterior for each merge.}
	\item{clust.number}{The number of clusters in each merge.}
	\item{cut}{The value of height corresponding to the maximum of the log posterior.}
	\item{transformed.par}{The transformed values of the model hyperparameters. 
	  The log transofmration is applied for variances, the idendtity for the mean, and the logit for the proportions.}
	\item{labels}{The labels associated to each clustering individual.}
	\item{effect.family}{The distribution assigned to the disappearing random effect in the function arguments.}
	\item{var.select}{The varible selection chosen in the function arguments.}
      }
     
\seealso{
		\link{loglikelihood}, \link{meancss}
	      }

	      
\examples{
data(gaelle)

# unreplicated clustering
gaelle.bclust<-bclust(x=gaelle,transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68)) 
par(mfrow=c(2,1))
plot(as.dendrogram(gaelle.bclust))
abline(h=gaelle.bclust$cut)
plot(gaelle.bclust$clust.number,gaelle.bclust$logposterior,xlab="",ylab="",type="b")
abline(h=max(gaelle.bclust$logposterior),xlab="",ylab="")

#replicated clustering
gaelle.id<-rep(1:14,c(3,rep(4,13))) # first 3 rows replication of ColWT , the other mutants each 
gaelle.lab<-c("ColWT","d172","d263","isa2",
"sex4","dpe2","mex1","sex3","pgm","sex1","WsWT","tpt","RLDWT","ke103")

gaelle.bclust<-bclust(gaelle,rep.id=gaelle.id,labels=gaelle.lab,transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
plot(as.dendrogram(gaelle.bclust))
abline(h=gaelle.bclust$cut)
plot(gaelle.bclust$clust.number,gaelle.bclust$logposterior,xlab="",ylab="",type="b")
abline(h=max(gaelle.bclust$logposterior),xlab="",ylab="")
}


