\name{GCM_runs}
\alias{GCM_runs}

\title{GCM runs under several forcing scenarios}
\description{outputs from general circulation models (GCMs) under several forcing scenarios}

\format{a list of two elements: the first element is a 2107 by 3 matrix with rows for the number of grid cells and columns for the number of model runs under anthropogenic forcing; the second element is a 2107 by 4 matrix with rows for the number of grid cells and columns for the number of model runs under natural forcing}

\examples{
# plot GCM runs for anthropogenic forcing and natural forcing scenarios
boxplot(GCM_runs[[1]], xlab="number of runs under anthropogenic forcing",
        ylab="model runs", cex=0.5, cex.axis=0.6, cex.lab=0.8, 
        pars=list(boxwex=0.5), medcol="red", medlwd=1)
boxplot(GCM_runs[[2]], xlab="number of runs under natural forcing",
        ylab="model runs", cex=0.5, cex.axis=0.6, cex.lab=0.8, 
        pars=list(boxwex=0.5), medcol="red", medlwd=1)
}
\keyword{datasets}
