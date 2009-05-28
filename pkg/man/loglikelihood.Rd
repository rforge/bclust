\name{loglikelihood}

\alias{loglikelihood}

\title{computes the model log likelihood}

\description{The function is useful for deriving the maximum likelihood
  estimates of the model parameters.}

\usage{loglikelihood(x.mean,x.css,repno,transformed.par,effect.family="gaussian",var.select=TRUE)}

 \arguments{
   \item{x.mean}{The mean matrix of the clustering individuals.}
   
   \item{x.css}{The corrected sum of squares matrix of the clustering individuals.}
   \item{repno}{The vector containing the number of replications of each
     clustering individual.}
   
 \item{transformed.par}{The vector of transformed model parameters in which the data likelihood will be evaluated at. The tranformation is the log for variance parameters, the identity for the mean, and the logit for the proportions. The length of the vector depends 
   on the chosen \code{effect.family} and \code{var.select}.}
 
 \item{effect.family}{Distribution family of the disappearing random
 components. Choices are "gaussian" or "alaplace" allowing Gaussian or
 asymmetric Laplace family, respectively. }

 \item{var.select}{A logcial value, TRUE for fitting models that define
   spike-and-slab in variable level.}
 
  	}

\details{Sometimes estimation of the model parameters is difficult,
always check the convergence of the optimisation algorithm. The
asymmetric Laplace model, \code{effect.family="alaplace"} is often more
difficult to optimise than \code{effect.family="gaussian"}.

If your data is standardised (having general mean
equal to zero and general variance equal one) often the
\code{loglikelihood} function is maximised over values between -5 and 5.

The \code{transformed.par} is a vector of transformed model parameters
having length 5 up to 7 depending on the chosen model.

The \code{transformed.par} is \code{(log sigma^2,log sigma_eta^2, log
sigma_theta^2,mu, logit p,logit q )} a vector of length 6 for
\code{effect.family="gaussian"} and \code{var.select=TRUE}, and is (log
sigma^2,log sigma_eta^2,log sigma_theta_L^2,log sigma_theta_R^2,mu,
logit p, logitq) a vector of length 7 for
\code{effect.family="alaplace"}  and \code{var.select=TRUE}.

When \code{var.select=FALSE} the \code{q} parameter is dropped yielding a vector
of length 5 for \code{effect.family="gaussian"} and a vector of length 6
for \code{effect.family="alaplace"}.

The model parameters has natural interpretations, \code{sigma^2} is the
variance of measurement error; \code{sigma^2_eta} is the variance of
between-individual error; \code{sigma^2_theta} is the variance of
disappearing random component which is decomposed to
\code{sigma^2_theta_L, sigma^2_theta_R} the left and the right tail
variances if the model is asymmetric Laplace; \code{mu} is the general
mean; \code{p} is the proportion of active variable-clusters, and
\code{q} is the proportion of active variables. For more details about
the models and parameters, see the website of the package:
\url{http://bclust.probstat.ch}. }

\seealso{
	\link{bclust},\link{meancss} 
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

