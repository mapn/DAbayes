#' @title Parallel inference in Bayesian hierarchical model for climate-change detection 
#' and attribution problem
#'
#' @param ensemble an n by N matrix of ensembles of measured variable, say temperature increase
#' 
#' @param model_runs a list of elements, each of which is model-output of the measured variable 
#'        (temperature increase) under a specific forcing scenario
#'
#' @param control_runs an n by L0 matrix, where n is the number of grid cells, 
#'        and L0 is the number of control runs
#'        
#' @param r_max an integer representing the maximum number of empirical orthogonal
#'        functions
#' 
#' @param theta_intVal a list of 3 elements containing initial values for parameters
#'        in the MCMC algorithm
#'
#' @param prior a list of 6 elements in prior distributions for corresponding
#'        parameters
#'
#' @param AD_proposal a list of 6 elements containing quantities to adjust mean and 
#'        covariance in the proposal distribution
#'
#' @param niter an integer containing the total number of MCMC interations
#' 
#' @param start_adapting an integer specifiying when to adapt proposal in 
#'        the MCMC algorithm
#'
#' @param numCore an integer specifying how many cores are used in parallel computing
#'
#' @return a list of r elements, where r is the number of EOFs specified, and
#' each element is again a list of 6 elements containing posterior quantities of 
#' parameters, log-likelihood, chisq statistics, and prior:
#' 
#' beta: a matrix holds the posterior samples for the parameter beta with each 
#'              row corresponding to each beta
#'              
#' logsigma: a vector holds the posterior samples for the parameter log of sigma
#' 
#' lambda: a vector holds the posterior samples for the parameter lambda         
#' 
#' loglik: a vector holds the log-likelihood evaluated with updated parameters
#' 
#' chisq: a vector holds chisquare statistics for residual consistency test
#' 
#' prior: a vector holds the prior density evaluated with updated parameters
#'          
#' @description This function wraps all the necessary functions in DAbayes package
#' to run the MCMC Suite. By default, this will execute the adaptive MCM algorithm
#' in single core (sequentially), and it can be executed in parallel by changing
#' the number of cores to the number of cores available.
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export 
#' 
#' @keywords models
#' 
#' 
#' @examples 
#' \dontrun{ 
#' ## Run Gibbs sampling with Metropolis-Hastings algorithm in the 
#' # Bayesian statistical model
#' data(ensemble_temperature)
#' data(GCM_runs)
#' data(GCM_control_run)
#' 
#' Res <- DAbayesSuite(ensemble=ensemble_temperature, model_runs=GCM_runs, 
#' control_runs=GCM_control_run, r_max=2, theta_intVal=NULL, prior=NULL,
#' AD_proposal=NULL, niter=20000, start_adapting=100, numCore=1)
#' 
#' ## diagnostic plots for posterior samples, and more sophisticated
#' # diagnostic plots can be obtained with coda package
#' # take the first element in Res for example
#' # plot only first 2 beta's
#' plot(Res[[1]]$beta[1, ], xlab="Iteration", ylab="beta_1",  type="l")
#' plot(Res[[1]]$beta[2, ], xlab="Iteration", ylab="beta_2", type="l")
#' plot(Res[[1]]$logsigma, xlab="Iteration", ylab="logsigma", type="l")
#' # plot only first 3 lambda's
#' plot(Res[[1]]$lambda[1, ], xlab="Iteration", ylab="lambda1", type="l")
#' plot(Res[[1]]$lambda[2, ], xlab="Iteration", ylab="lambda2", type="l")
#' plot(Res[[1]]$lambda[3, ], xlab="Iteration", ylab="lambda3", type="l")
#' 
#' ## save plots 
#' burnin <- 5000
#' plotMCMC(Res, N, n, dir=getwd(), burnin)
#' }
#' 




###########################################################
###########################################################
DAbayesSuite <- function(ensemble, model_runs, control_runs, r_max=NULL,
           theta_intVal=NULL, prior=NULL, AD_proposal=NULL, niter=20000, 
           start_adapting=50, numCore=1){

### read data
y <- ensemble
x <- model_runs
x0 <- control_runs

###########################################################
###########################################################
### read dimensions of input data
m <- length(x) # number of forcing scenarios
n <- dim(y)[1] # number of spatial grid cells on the globe
N <- dim(y)[2] # number of temperature ensemble members
Lj <- sapply(x, ncol) # number of runs for each scenario

# Check that there are enough control runs
if(is.null(r_max)){
  r_max <- dim(control_runs)[2]-1
}

if(r_max>dim(x0)[2]-1){
  # include a warning message
  r_max <- dim(x0)[2]-1  
  warning("The range of EOFs exceeds the number of control runs!")
}

r_all <- 1:r_max  # specify range of EOFs

# estimate PCs and eigenvalues
empiricalcovmat <- cov(t(x0))

# compute first r eigenvalues and eigenvectors
temp <- RSpectra::eigs_sym(empiricalcovmat, max(r_all), which="LM")
B_all <- temp$vectors
K_all <- temp$values

# center the data
y <- y - mean(y)
for(j in 1:m){
  x[[j]] <- x[[j]]-mean(x[[j]])
}

### precomputation for W
outW <- parcomputeW(y,x)





###############################################################
###############################################################
###############################################################



### define quantities for BMA
r_num <- length(r_all)
loglik <- matrix(0, niter, r_num)
chisq_statistic <- matrix(0, niter, r_num)
prior_all <- matrix(0, niter, r_num)




### MCMC output quantities and initial values of parameters



### quantities in prior distribution
if (is.null(prior)){
  beta_mean <- matrix(0,m,1)
  beta_var <- diag(rep(1,m))
  divisionfactor_lambda <- 3^2
  lambda_mean_all <- log(K_all)
  lambda_var_all <- var(lambda_mean_all)/divisionfactor_lambda*diag(rep(1,max(r_all)))
  logsigma_mean <- log((sum(diag(empiricalcovmat)) - sum(lambda_mean_all))/n)
  logsigma_var <- 1
}else{
  beta_mean <- prior$beta_mean
  beta_var <- prior$beta_var
  lambda_mean_all <- prior$lambda_mean_all
  lambda_var_all <- prior$lambda_var_all
  logsigma_mean <- prior$logsigma_mean
  logsigma_var <- prior$logsigma_var
}


### initial values for paramters
if (is.null(theta_intVal)){
  beta <- rep(0,m)
  lambda_all <- vector("list", r_num)
  lambda_initial <- vector("list", r_num)
  for(r_ind in 1:r_num){
    lambda_initial[[r_ind]] <- lambda_mean_all[1:r_ind]
  }
  logsigma <- logsigma_mean
}else{
  beta <- theta_intVal$beta
  lambda_initial <- theta_intVal$lambda
  logsigma <- theta_intVal$logsigma
}


### initial values for proposal distributions
if(is.null(AD_proposal)){
  multiplier_cov_beta_temp <- 0.01
  mean_beta <- matrix(0,m,1)
  cov_beta <- multiplier_cov_beta_temp^2*diag(rep(1,m))
  mean_logsigma <-logsigma[1]
  cov_logsigma <- logsigma_var
}else{
  mean_beta <- AD_proposal$mean_beta
  cov_beta <- AD_proposal$cov_beta
  mean_logsigma <- AD_proposal$mean_logsigma
  cov_logsigma <- AD_proposal$cov_logsigma
}

### begin parallel for-loop

# register cores
doMC::registerDoMC(numCore)
# check number of registered cores
if (!getDoParWorkers()){
  sprintf("Cores are registered successfully", )
}


# define initial values
theta <- vector("list", 3)
names(theta) <- c("beta","logsigma","lambda")
theta$beta <- beta
theta$logsigma <- logsigma

prior_theta <- vector("list", 6)
names(prior_theta) <- c("beta_mean","beta_var","lambda_mean","lambda_var",
                          "logsigma_mean","logsigma_var")
prior_theta$beta_mean <- beta_mean
prior_theta$beta_var <- beta_var
prior_theta$logsigma_mean <- logsigma_mean
prior_theta$logsigma_var <- logsigma_var

ADtheta <- vector("list", 6)
names(ADtheta) <- c("mean_beta","cov_beta","mean_logsigma",
                      "cov_logsigma","mean_lambda","cov_lambda")
ADtheta$mean_beta <- mean_beta
ADtheta$cov_beta <- cov_beta
ADtheta$mean_logsigma <- mean_logsigma
ADtheta$cov_logsigma <- cov_logsigma
      

Results <- foreach(r_ind=1:r_num) %dorng% {
  
  
  # quantities that depend on r explicitely
  r <- r_all[r_ind]
  B <- as.matrix(B_all[ ,1:r],n,r)
  
  outW$B <- B
  
  lambda_mean <- lambda_mean_all[1:r]
  
  
  # initial values for parameters in MCMC
  lambda <- as.matrix(lambda_initial[[r_ind]], r_ind, 1)
  theta$lambda <- lambda


  # quantities in prior distribution
  lambda_mean <- as.matrix(lambda_mean_all[1:r_ind], r_ind, 1)
  lambda_var <- lambda_var_all[1:r_ind,1:r_ind]
  prior_theta$lambda_mean <- lambda_mean
  prior_theta$lambda_var <- lambda_var



  # quantities in proposal distributions
  mean_lambda <- lambda_initial[[r_ind]]
  cov_lambda <- lambda_var
  ADtheta$mean_lambda <- mean_lambda
  ADtheta$cov_lambda <- cov_lambda


  # compute likelihood 
  out.loglik <- parDA_loglik(outW,theta,N,Lj)
  loglik_old <- out.loglik$loglik
  chisq_old <- out.loglik$chisq
  
  Res <- parDA_GIBBS(outW,theta,prior_theta,ADtheta,N,Lj,niter,start_adapting)
  
  return(Res)
}

return(Results)

}


