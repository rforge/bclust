\name{loglikelihood}

\alias{loglikelihood}

\title{computes the model log likelihood useful for estimation of the transformed.par}

\description{The function is useful for deriving the maximum likelihood
  estimates of the model parameters.}

\usage{loglikelihood(x.mean,x.css,repno,transformed.par,
effect.family="gaussian",var.select=TRUE)}

 \arguments{
   \item{x.mean}{The mean matrix of the clustering types from the \link{meancss} function.}
   
   \item{x.css}{The corrected sum of squares matrix of the clustering types from the \link{meancss} function.}
   \item{repno}{The vector containing the number of replications of each
     clustering type corresponding to the each row of \code{x.mean} and \code{x.css}, from the \link{meancss} function.}
   
 \item{transformed.par}{The vector of transformed model parameters that the data likelihood will be evaluated at. The transformation is the log for the variance parameters, the identity for the mean, 
 and the logit for the proportions. The length of the vector depends 
   on the chosen \code{effect.family} and \code{var.select}.}
 
 \item{effect.family}{Distribution family of the disappearing random
 components. Choices are "gaussian" or "alaplace" allowing Gaussian or
 asymmetric Laplace family, respectively. }

 \item{var.select}{A logical value, \code{TRUE} for fitting models that define
   spike-and-slab in variable level, thus allowing Bayesian variable selection.}
   	}

\details{Sometimes estimation of the model parameters is difficult,
always check the convergence of the optimisation algorithm. The
asymmetric Laplace model, \code{effect.family="alaplace"}, is often more
difficult to optimise than \code{effect.family="gaussian"}.

If your data is standardised (having general mean zero and general variance one) the
\code{loglikelihood} function is usually maximised over values between -5 and 5.

The \code{transformed.par} is a vector of transformed model parameters
having length 5 up to 7 depending on the chosen model.

The \code{transformed.par} is \eqn{(\log \sigma^2,\log \sigma_\eta^2, \log
\sigma_\theta^2,\mu, \mathrm{logit} p,\mathrm{logit} q }{log s2, log s2_h, log s2_t, m, logit p, logit q} a vector of length 6 for \code{effect.family="gaussian"} and \code{var.select=TRUE}, and is \eqn{(\log
\sigma^2,\log \sigma_\eta^2,\log \sigma_{\theta_L}^2,\log \sigma_{\theta_R}^2,\mu,
\mathrm{logit} p, \mathrm{logit} q)}{log s2, log s2_h, log s2_tL, log s2_tR, m, logit p, logit q} a vector of length 7 for \code{effect.family="alaplace"}  and \code{var.select=TRUE}.

When \code{var.select=FALSE} the \eqn{q} parameter is dropped, yielding a vector
of length 5 for \code{effect.family="gaussian"} and a vector of length 6
for \code{effect.family="alaplace"}.

We assumed a Bayesian linear model being \deqn{y_{vctr}=\mu+\eta_{vct}+\delta_v \gamma_{vc}\theta_{vc}+\varepsilon_{vctr}}{y_{vctr}=m+h_{vct}+d_{v}*g_{vc}*t_{vc}+e_{vctr}} where \eqn{y_{vctr}} is the available data on variable \eqn{v}, 
cluster(or class) \eqn{c}, type \eqn{t}, and replicate \eqn{r}; \eqn{\eta_{vct}}{h_{vct}}
is the between-type error, \eqn{\theta_{vc}}{t_{vc}} is the disappearing random component controlled by the Bernoulli variables \eqn{\delta_v}{d_{v}} with success probability \eqn{q} and \eqn{\gamma_{vc}}{g_{vc}} with 
success probability \eqn{p}; and \eqn{\varepsilon_{vctr}}{e_{vctr}} is the between-replicate error. The types inside a cluster (or class) share the same \eqn{\theta_{vc}}{t_{vc}}, but may arise with a different  \eqn{\eta_{vct}}{h_{vct}}, 
for more details see the package website \url{http://bclust.probstat.ch} and documents there in.

The model parameters has natural interpretations, \eqn{\sigma^2}{s2} is the
between replicate error variance; \eqn{\sigma^2_\eta}{s2_h} is the variance of
between-type error; \eqn{\sigma^2_\theta}{s2_t} is the variance of
the disappearing random component which is decomposed to
\eqn{\sigma^2_{\theta_L}}{s2_tL}, \eqn{\sigma^2_{\theta_R}}{s2_tR} the left and the right tail
variances if the model is asymmetric Laplace; \eqn{\mu}{m} is the general
level; \eqn{p} is the proportion of active variable-cluster (or variable-class) combinations, and
\eqn{q} is the proportion of the active variables. For more details about
the model and parameters, see the website of the package \url{http://bclust.probstat.ch}. }

\seealso{
	\link{bclust},\link{bdiscrim}, \link{meancss}.
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

transpar<-optim(rep(0,6),optimfunc,method="BFGS")$par #gives argmin(-loglikelihood)
#Attention: put a vector of correct length for the evaluation of the likelihood

plot(bclust(gaelle,transformed.par=transpar))

        }

