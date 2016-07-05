context("testing function computeW()")

test_that("computeW", {

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

neighborhood_matrix <- Matrix::sparseMatrix(1:n, 1:n, x=rep(1,n))


outW <- computeW(y,x,B,neighborhood_matrix,gamma_prior)


expect_equal(length(outW), 7)

expect_equal(length(outW[[1]]), ng)
expect_equal(length(outW[[2]]), ng)
expect_equal(is.list(outW[[3]]), TRUE)
expect_equal(length(outW[[3]]), ng)
expect_equal(is.list(outW[[4]]), TRUE)
expect_equal(length(outW[[4]]), ng)
expect_equal(identical(dim(outW[[5]]), as.integer(c(m,ng))), TRUE)
expect_equal(identical(dim(outW[[6]]), as.integer(c(r,ng))), TRUE)
expect_equal(is.list(outW[[7]]), TRUE)
expect_equal(length(outW[[7]]), ng)


})
