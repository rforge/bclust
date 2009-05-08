 \name{bclust}
\alias{bclust}
\title{Bayesian agglomerative clustering for high dimensional data with variable selection.}
\description{The function clusters data saved in a matrix using an additive linear model with disappearing random effects. 
The model has built-in spike-and-slab components which quantifies important variables for clustering.}

(x,model.par,x.id=1:nrow(x),effect.family="gaussian",
var.select=TRUE,labels=NULL)

\usage{bclust(x,model.par,x.id=1:nrow(x),effect.family="gaussian",
var.select=TRUE,labels=NULL)	}
 \arguments{
 \item{x}{A numeric matrix of clustering indiviuals in rows and variables in columns.}
 \item{model.par}{Model parameters in a vector. The length of the vector depends on the chosen model and the availability of variable selection.}
 \item{x.id}{A vector consisting of positive integer elements having the same length as the number of rows of \code{x}, identifying replicates of a type.}
 \item{effect.family}{Distribution family of the disappearing random components. Choices are "gaussian" or "alaplace" allowing Gaussian or asymmetric Laplace family, respectively. } 
 \item{var.select}{A logcial value, TRUE for fitting models that define spike-and-slab in variable level.}
 \item{labels}{A vector of strings refering to the labels of clustering individuals. The length of the vector should match to \code{max(x.id)}}
 	}
\details{The function calls internal C functions depending on the chosen model. The C-stack  of the system may overflow if you have a large dataset. You may need to adjust the stack before running 
R using your operational system command line. 
		}
\value{
	\item{data}{The data matrix reordered according to \code{x.id}}
	\item{repno}{the number of replicates of the values of \code{x.id}}
	\item{merge}{The merge matrix \code{hclust} format.}
	\item{height}{A monotone vector reffering to height of dendrogram.}
	\item{logposterior}{The log posterior for each merge.}
	\item{clust.number}{The number of clusters for each merge.}
	\item{cut}{The value of height corresponding to the maximum of the log posterior.}
	\item{transformed.par}{The values of \code{model.par} after transformation. The log transofmration is applied for variances, the idendtity for the mean, and the logit for the proportions.}
	\item{labels}{The labels associated to each }
	}
\seealso{
		\link{likelihood}, \link{meancss}
			}
\examples{
data(gaelle)
bclust.gaelle<-bclust(x=gaelle,transformed.par=c(-1.8,-1.7,0.2,0.2,0.6,0)) #transformed.par chosen close to the values obtained by maximum likelihood
par(mfrow=c(2,1))
plot(bclust.gaelle)
abline(h=bclust.gaelle$cut)
plot(bclust.gaelle$clust.number,bclust.gaelle$logposterior,xlab="",ylab="",type="b")
abline(h=max(bclust.gaelle$logposterior),xlab="",ylab="")
}

