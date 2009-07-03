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
     clustering individual corresponding to each row of \code{x.mean} and \code{x.css}.}
   
 \item{transformed.par}{The vector of transformed model parameters that the data likelihood will be evaluated at. The tranformation is the log for variance parameters, the identity for the mean, 
 and the logit for the proportions. The length of the vector depends 
   on the chosen \code{effect.family} and \code{var.select}.}
 
 \item{effect.family}{Distribution family of the disappearing random
 components. Choices are "gaussian" or "alaplace" allowing Gaussian or
 asymmetric Laplace family, respectively. }

 \item{var.select}{A logcial value, \code{TRUE} for fitting models that define
   spike-and-slab in variable level, thus allowing Bayesian variable selection.}
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

The \code{transformed.par} is \eqn{(\log \sigma^2,\log \sigma_\eta^2, \log
sigma\_\theta^2,\mu, logit p,logit q } a vector of length 6 for
\code{effect.family="gaussian"} and \code{var.select=TRUE}, and is \eqn{(\log
\sigma^2,\log \sigma_\eta^2,\log \sigma_{\theta_L}^2,\log \sigma_{\theta_R}^2,\mu,
logit p, logit q)} a vector of length 7 for
\code{effect.family="alaplace"}  and \code{var.select=TRUE}.

When \code{var.select=FALSE} the \code{q} parameter is dropped yielding a vector
of length 5 for \code{effect.family="gaussian"} and a vector of length 6
for \code{effect.family="alaplace"}.

The assumed Bayesian model for clustering is \deqn{y_{vctr}=\mu+\eta_{vct}+\delta_v \gamma_{vc}\theta_{vc}+\varepsilon_{vctr};} where \eqn{y} is the avilable data on variable \eqn{v}, 
initial cluster \eqn{c}, type \eqn{t}, and replicate \eqn{r}; \eqn{\eta} 
is between type error, \eqn{\theta} is the disppearing random component controlled by bernoulli variables \eqn{\delta} with succes probability \eqn{q} and \eqn{\gamma} with 
success probability \eqn{p}; and \eqn{\varepsilon} is the between replicate error.

The model parameters has natural interpretations, \eqn{\sigma^2} is the
between replicate error variance; \eqn{\sigma^2_\eta} is the variance of
between-type error; \eqn{\sigma^2_\theta} is the variance of
the disappearing random component which is decomposed to
\eqn{\sigma^2_{\theta_L}, \sigma^2_{\theta_R}} the left and the right tail
variances if the model is asymmetric Laplace; \eqn{\mu} is the general
mean; \eqn{p} is the proportion of active variable-cluster combinations, and
\eqn{q} is the proportion of the active variables. For more details about
the model and parameters, see the website of the package:
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

transpar<-optim(rep(0,6),optimfunc)$par #gives argmin(-loglikelihood)
#Attention should be paid to put a vector of correct length for evaluation of the likelihood

plot(bclust(gaelle,transformed.par=transpar))

        }

