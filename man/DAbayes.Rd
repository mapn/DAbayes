\name{DAbayes}
\alias{DAbayes}
\title{DAbayes: Climate-Change Detection and Attribution Based on the Bayesian Hierarchical Model}
\description{Climate change is happening, and human forcing is suspected to be one of the causes. By considering climate-change detection and attribution as a multivariate spatial or spatio-temporal regression problem, we would like to understand whether and how the climate has changed due to anthropogenic and natural forcing scenarios. A Bayesian hierarchical statistical model has been developed for the detection and attribution problem, and the DAbayes package implements the Bayesian statistical model, which allows sequential and parallel Bayesian inference (e.g., Markov Chain Monte Carlo (MCMC) algorithm) via Bayesian model averaging. This package incorporates both a simulation example and real data analysis. For the latter, observed or reconstructed measurements (e.g., temperature change) and general circulation model (GCM) outputs under several forcing scenarios are used to illustrate the statistical problem.}


\details{
The DAbayes package implements a Bayesian statistical model with adaptive MCMC algorithm run in sequential (with single core) or parallel (with multiple cores).

\bold{Major methods}

\code{\link{DAbayesSuite}} This function runs entire inference algorithm and returns the posterir samples, which can be used for further diagonostic and results reporting.

\code{\link{update_loglik}} and \code{\link{parDA_loglik}} These functions compute the loglikelihood function in the model in sequential and parallel algorithms.

\code{\link{DA_theta_MH}} and \code{\link{parDA_loglik}} These functions update parameters via Gibbs sampling with Metropolis-Hasting algortihm in sequential and parallel algorithms.


\code{\link{update_var}} This function adaptively adjusts parameters in proposal distributions.


\bold{Other functions}

\code{\link{DA_GIBBS}} and \code{\link{parDA_GIBBS}} These functions update parameters for many iterations and adaptively adjust parameters in proposal distributions, which are used in sequential and parallel algorithms.

\code{\link{DA_ADtheta_update}} This function adaptively adjusts parameters in proposal distributions for given updated parameters in the MCMC algorithm.

\bold{Datasets}

\code{\link{ensemble_temperature}} A dataset containing linear trend from an ensemble of observed or reconstructed temperatures.

\code{\link{GCM_runs}} outputs from general circulation models (GCMs) under several forcing scenarios.

\code{\link{GCM_control_run}} Outputs from general circulation models (GCMs) without any external forcing represent internal climate variability and GCM variability.



\bold{Examples}\preformatted{
## Run Gibbs sampling with Metropolis-Hastings algorithm in the 
# Bayesian statistical model
library(DAbayes)
data(ensemble_temperature)
data(GCM_runs)
data(GCM_control_run)

Res <- DAbayesSuite(ensemble=ensemble_temperature, model_runs=GCM_runs, 
control_runs=GCM_control_run, r_max=2, theta_intVal=NULL, prior=NULL,
AD_proposal=NULL, niter=5000, start_adapting=50, numCore=4)

## diagnostic plots for posterior samples, and more sophisticated
# diagnostic plots can be obtained with coda package
# take the first element in Res for example
# plot only first 2 beta's
plot(Res[[1]]$beta[1, ], xlab="Iteration", ylab="beta_1",  type="l")
plot(Res[[1]]$beta[2, ], xlab="Iteration", ylab="beta_2", type="l")
plot(Res[[1]]$logsigma, xlab="Iteration", ylab="logsigma", type="l")
# plot only first 3 lambda's
plot(Res[[1]]$lambda[1, ], xlab="Iteration", ylab="lambda1", type="l")
plot(Res[[1]]$lambda[2, ], xlab="Iteration", ylab="lambda2", type="l")
plot(Res[[1]]$lambda[3, ], xlab="Iteration", ylab="lambda3", type="l")

## save plots 
burnin <- 1000
plotMCMC(Res, N, n, dir=getwd(), burnin)

} % end preformatted
}

\concept{Detection and Attribution, ensemble, ensemble member, model output, general circulation model (GCM), Bayesian Statistical Model,
Markov Chain Monte Carlo Algorithm}

\note{
This package can run the Bayesian inference (MCMC algorithm), and save several crucial plots for D&A problem, but there is no function implemented yet to make diagnostics based on posterior samples, which can be obtained through coda package.
}
