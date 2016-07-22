\name{GCM_control_run}
\alias{GCM_control_run}
\title{GCM control run without any external forcing}

\description{Outputs from general circulation models (GCMs) without any external forcing represent internal climate variability and GCM variability}

\format{ an n by L0 matrix, where n is the number of grid cells, and L0 is the number of control runs.}

\examples{
# plot GCM control runs 
boxplot(GCM_control_run, xlab="number of runs without any external forcing scenario",
        ylab="control runs", cex=0.5, cex.axis=0.6, cex.lab=0.8, 
        pars=list(boxwex=0.5), medcol="red", medlwd=1)
}
\keyword{datasets}
