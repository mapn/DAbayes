#' @title Simulate data to demonstrate the Bayesian statistical model
#' 
#' @param n an integer containing the number of spatial grid cells over the globe
#' 
#' @param N an integer containing number of ensemble members for measured variable
#' 
#' @param m an integer containing the number of forcing scenarios
#' 
#' @param Lj a vector containing the number of model control runs under each forcing
#' scenarios
#' 
#' @param L0 an integer containing the number of model control runs without any
#'  external forcing
#' 
#' @param trend trend of the simulated data for the ensemble of the measured variable
#' with default value 30
#' 
#' @return a list of 3 elements containing an n by N matrix of the ensemble of the
#' measured variable, an n by m matrix of model outputs for measured variable under
#' m forcing scenarios, and an n by L0 matrix of model control runs for measured  
#' variable without any external forcing scenario
#' 
#' @description This function simulates the ensemble of the measured variable, model 
#' outputs for the measured variable under different forcing scenarios, and model outputs
#' for the measured variable without any external forcing scenario.
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @keywords models
#' 
#' @export
#' 
#' @examples 
#' n <- 30
#' N <- 10
#' m <- 3
#' Lj <- c(5, 3, 7)
#' L0 <- 8
#' trend <- 30
#' DAdata <- simDAdata(n, N, m, Lj, L0, trend)
#' # ensembles of the measured variable
#' y <- DAdata[[1]]
#' # model outputs for the measured variable under different forcing scenarios
#' x <- DAdata[[2]]
#' # model outputs for the measured variable without any external forcing scenario
#' x0 <- DAdata[[3]]



simDAdata <- function(n=30, N=10, m=3, Lj=c(5,3,7), L0=8, trend=30){
  
  # matrix containing ensembles of measured temperature 
  # rows correspond to different grid cells on the globe
  y <- matrix(rnorm(n*N, 0, 1), nrow=n, ncol=N) + trend
  
  # number of ensemble runs for each forcing scenario
  numforcings_sim <- Lj;
  
  # model-output temperature increases for each forcing scenario j
  x <- vector(mode = "list", 3);
  for(j in 1:length(x))
  {
    x[[j]] <- matrix(rnorm(n*numforcings_sim[j]), 
                     nrow=n, ncol=numforcings_sim[j])
  }
  
  # temperature increases in control runs
  x0 <- matrix(rnorm(n*L0), nrow=n, ncol=L0)
  
  return(list(y, x, x0))
}
  
