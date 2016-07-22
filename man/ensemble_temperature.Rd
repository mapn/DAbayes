\name{ensemble_temperature}
\alias{ensemble_temperature}
\title{Observed or reconstructed temperatures (e.g., from CRU or NCDC)}

\description{A dataset containing linear trend from an ensemble of observed or reconstructed temperatures} 

\format{ A n by N matrix, where n=2107 is the number of grid cells over globe, and N=394 is the number of ensemble memembers}

\examples{
# estimate linear trend
\dontrun{
data("ensemble_temperature")
data("GCM_runs")
y <- ensemble_temperature
x <- GCM_runs
N <- dim(y)[2]
avg <- function(x) {return(apply(x, 1, mean))}
xbar <- sapply(x, avg)
Y <- c(y)
X <- kronecker(diag(rep(1,N),N,1), xbar)

# OLS estimates
fit <- lm(Y~X)
summary(fit)
}

# plot the first 10 ensemble memebers
boxplot(ensemble_temperature[ ,1:10], xlab="ensemble member", 
        ylab="temperature increase", cex=0.5, cex.axis=0.6, cex.lab=0.8, 
        pars=list(boxwex=0.5), medcol="red", medlwd=1)
}
\keyword{datasets}
