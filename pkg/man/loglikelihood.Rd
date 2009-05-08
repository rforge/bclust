 \name{loglikelihood}
\alias{loglikelihood}
\title{computes the model log likelihood}
\description{The function is useful for deriving the maximum likelihood estimates of the model parameters.}
\usage{loglikelihood(x.mean,x.css,repno,transformed.par,effect.family="gaussian",var.select=TRUE)}
 \arguments{
 \item{x.mean}{The mean matrix of the clustering individuals.}
 \item{x.css}{The corrected sum of squares matrix of the clustering individuals.}
 \item{repno}{The vector containing the number of replications of each clustering individual.}
 \item{transformed.par}{The vector of transformed model parameters in which the data likelihood will be evaluated at. The tranformation is the log for variance parameters, the identity for the mean, and the logit for the proportions. The length of the vector depends 
 on the chosen \code{effect.family} and \code{var.select}.}
 \item{effect.family}{Distribution family of the disappearing random components. Choices are "gaussian" or "alaplace" allowing Gaussian or asymmetric Laplace family, respectively. }
 \item{var.select}{A logcial value, TRUE for fitting models that define spike-and-slab in variable level.}
  	}
\details{Sometimes estimation of the model parameters is difficult, always check the convergence of the optimisation algorithm. The asymmetric Laplace model, \code{effect.family="alaplace"} is often more difficult
to optimise. If your data is standardised (having general mean equal to zero and general variance equal one) often the \code{loglikelihood} function is 
maximised over values between -5 and 5.
The \code{transformed.par} is a vector of transformed model parameters having length 5 up to 7 depending on the chosen model.
The \code{transformed.par} is $(\log \sigma^2,\log \sigma_\eta^2,\log sigma_\theta^2,\mu,\mathrm{logit}p,\mathrm{logit}q)$ a vector of length 6 for \code{effect.family="gaussian"} 
and \code{var.select=TRUE}, and is $(\log \sigma^2,\log \sigma_\eta^2,\log sigma_{\theta_L}^2,\log sigma_{\theta_R}^2,\mu,\mathrm{logit}p,\mathrm{logit}q)$ 
a vector of length 7 for \code{effect.family="alaplace"}. When \code{var.select=FALSE} the parameter is dropped yielding a vector of length 5 for \code{effect.family="gaussian"}  and a vector of 
length 6 for \code{effect.family="alaplace"}. The model parameters has natural interpretations, $\sigma^2$ is the variance of measurement error, 
$\sigma^2_\eta$ is the variance of between-individual error, $\sigma^2\theta$ is the variance of disappearing random component 
which is decomposed to $\sigma^2_{\theta_L}, \sigma^2_{\theta_R}$  the left and the right tail variances if the asymmetric Laplace. $\mu$ is the general mean, $p$ is the proportion of active variable-clusters, and
$q$ is the proportion of active variables. For more details see the website of the package. }
\seealso{
		\link{meancss}, \link{bclust}, 
			}
\examples{
data(gaelle)
gaelle.id<-rep(1:14,c(3,rep(4,13))) # first 3 rows replication of ColWT , the other mutants each have four replications
mc.gaelle<-meancss(gaelle,gaelle.id)

optimfunc<-function(theta)
{
-loglikelihood(x.mean=mc.gaelle$mean,x.css=mc.gaelle$css,
repno=mc.gaelle$repno,transformed.par=theta)#compute - log likelihood
}

transpar<-optim(rep(0,6),optimfunc)$par #gives argmin (- log like)
plot(bclust(gaelle,transformed.par=transpar))
           }

