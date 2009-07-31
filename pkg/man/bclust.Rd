\name{bclust}

\alias{bclust}

\title{Bayesian agglomerative clustering for high dimensional data with variable selection.}

\description{The function clusters data saved in a matrix using an additive linear model with disappearing random effects. 
The model has built-in spike-and-slab components which quantifies important variables for clustering and can be extracted using the \code{imp} function.}

\usage{bclust(x,rep.id=1:nrow(x),effect.family="gaussian",
var.select=TRUE,transformed.par,labels=NULL)
}


\arguments{
  \item{x}{A numeric matrix, with clustering individuals in rows
    and the variables in columns.}
 
\item{rep.id}{A vector consisting of positive integer elements having the
 same length as the number of rows of \code{x}. This vector identifies replicates of a clustering type such that the total number of clustering types is \code{max(rep.id)}. If nothing is declared the function
 presupposes that the data are unreplicated, that is each row of \code{x} is a clustering type.}
\item{effect.family}{Distribution family of the disappearing random components. Choices are "gaussian" or "alaplace"
 allowing Gaussian or asymmetric Laplace family, respectively. }

 \item{var.select}{A logical value, \code{TRUE} for fitting models that define
   spike-and-slab distribution in variable level and allows Bayesian variable selection.}

  \item{transformed.par}{The transformed
 model parameters in a vector. The length of the vector depends on the
 chosen model and the availability of variable selection. The log transformation is supposed to be applied for the variance parameters, the identity for the mean, and the logit for the proportions. The function \link{loglikelihood} can be used to estimate them from the data.}

\item{labels}{A vector of strings
 referring to the labels of clustering types. The length of the
 vector should match to \code{max(rep.id)}. The first element corresponds to the label of the type having the smallest integer value in \code{rep.id}, the second element refers to 
 the label of the type having the second smallest integer in \code{rep.id}, and so on.}} 

\details{The function calls internal \code{C} functions depending on the chosen
model. The C-stack of the system may overflow if you have a large
dataset. You may need to adjust the stack before running \code{R} using your
operation system command line. If you use Linux, open a console and type
\code{>ulimit -s unlimited}, then run \code{R} in the same console. The Microsoft Windows users don't need to increase the stack size.

We assumed a Bayesian linear model for clustering being \deqn{y_{vctr}=\mu+\eta_{vct}+\delta_v \gamma_{vc}\theta_{vc}+\varepsilon_{vctr}}{y_{vctr}=m+h_{vct}+d_{v}*g_{vc}*t_{vc}+e_{vctr}} where \eqn{y_{vctr}} is the available data on variable \eqn{v}, 
cluster \eqn{c}, clustering type \eqn{t}, and replicate \eqn{r}; \eqn{\eta_{vct}}{h_{vct}}
is the between-type error, \eqn{\theta_{vc}}{t_{vc}} is the disappearing random component controlled by the Bernoulli variables \eqn{\delta_v}{d_{v}} with success probability \eqn{q} and \eqn{\gamma_{vc}}{g_{vc}} with 
success probability \eqn{p}; and \eqn{\varepsilon_{vctr}}{e_{vctr}} is the between-replicate error. The types inside a cluster share the same \eqn{\theta_{vc}}{t_{vc}}, but may arise with a different  \eqn{\eta_{vct}}{h_{vct}}, 
for more details see the package website \url{http://bclust.probstat.ch} and documents there in.}

\value{
	\item{data}{The data matrix, reordered according to \code{rep.id}.}
	\item{repno}{The number of replicates of the values of \code{rep.id}}
	\item{merge}{The merge matrix, in \code{hclust} object format.}
	\item{height}{A monotone vector referring to the height of the
	constructed tree.}
	\item{logposterior}{The log posterior for each merge.}
	\item{clust.number}{The number of clusters for each merge.}
	\item{cut}{The value of the height corresponding to the maximum of the log posterior in agglomerative path.}
	\item{transformed.par}{The transformed values of the model parameters. 
	  The log transformation is applied for the variance parameters, the identity for the mean, and the logit for the proportions.}
	\item{labels}{The labels associated to each clustering type.}
	\item{effect.family}{The distribution assigned to the disappearing random effect in the function arguments.}
	\item{var.select}{The variable selection chosen in the function arguments.}
      }
     
\seealso{
		\link{loglikelihood}, \link{meancss}, \link{imp}.
	      }

	      
\examples{
data(gaelle)

# unreplicated clustering
gaelle.bclust<-bclust(x=gaelle,transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68)) 
par(mfrow=c(2,1))
plot(as.dendrogram(gaelle.bclust))
abline(h=gaelle.bclust$cut)
plot(gaelle.bclust$clust.number,gaelle.bclust$logposterior,xlab="Number of clusters",ylab="Log posterior",type="b")
abline(h=max(gaelle.bclust$logposterior))

#replicated clustering
gaelle.id<-rep(1:14,c(3,rep(4,13))) # first 3 rows replication of ColWT , the other mutants each 
gaelle.lab<-c("ColWT","d172","d263","isa2",
"sex4","dpe2","mex1","sex3","pgm","sex1","WsWT","tpt","RLDWT","ke103")

gaelle.bclust<-bclust(gaelle,rep.id=gaelle.id,labels=gaelle.lab,transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
plot(as.dendrogram(gaelle.bclust))
abline(h=gaelle.bclust$cut)
plot(gaelle.bclust$clust.number,gaelle.bclust$logposterior,xlab="Number of clusters",ylab="Log posterior",type="b")
abline(h=max(gaelle.bclust$logposterior))
}


