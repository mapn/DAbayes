context("testing function update_loglik()")

test_that("update_loglik",{
	
set.seed(1234)
# matrix containing ensemblee of measured temperature increases
# rows correspond to different grid cells on the globe
y <- matrix(rnorm(30*10, 0, 1), nrow=30, ncol=10)

# number of ensemble runs for each forcing scenario
numforcings_sim <- c(5,3,7);

# model-output temperature increases for each forcing scenario j
x <- vector(mode = "list", 3);
for(j in 1:length(x))
{
  x[[j]] <- matrix(rnorm(nrow(y)*numforcings_sim[j]), 
                   nrow=nrow(y), ncol=numforcings_sim[j])
}

# temperature increases in control runs
x0 <- matrix(rnorm(nrow(y)*8), nrow=nrow(y), ncol=8)

#################### end of simulation ####################


###########################################################
### read dimensions of input data
m <- length(x) # number of forcing scenarios
n <- dim(y)[1] # number of spatial grid cells on the globe
N <- dim(y)[2] # number of temperature ensemble members
Lj <- sapply(x, ncol) # number of runs for each scenario

### preprocessing 
# construct basis function matrix B with principal components
r <- 3
empiricalcovmat <- cov(t(x0))

# load rARPACK library to compute first r eigenvalues and eigenvectors
#library(rARPACK)
temp <- rARPACK::eigs_sym(empiricalcovmat, r, which="LM")
B <- temp$vectors
K_hat <- temp$values


lambda_mean <- log(K_hat)
lambda_var <- var(lambda_mean)/3^2


### pre-computation for W

# possible values for the discrete prior of gamma
gamma_prior <- c(0.5,0.99,1,1.01,2,5)
ng <- length(gamma_prior)

# neighborhood matrix (replace by real one)
neighborhood_matrix <- Matrix::sparseMatrix(1:n, 1:n, x=rep(1,n))

# obtain outW from function computeW # first parameter in update_loglik
outW <- computeW(y,x,B,neighborhood_matrix,gamma_prior)

# prepare theta # second parameter in update_loglik
beta <- rep(0,m)
lambda <- rep(0,r)
logsigma <- 0
gamma_ind <- 0

theta <- list(gamma_ind,beta,logsigma,lambda)
names(theta) <- c("gind","beta","logsigma","lambda")

loglik_gamma <- rep(0, ng)

for(j in 1:ng)
{
    theta$gind <- j
	loglik_gamma[j] <- update_loglik(outW,theta,n,N,Lj)
}

## begin testing
expect_equal(length(loglik_gamma), ng)
expect_equal(is.vector(loglik_gamma), TRUE)

})
