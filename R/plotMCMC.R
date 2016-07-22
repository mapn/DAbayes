#' @title Plot results from posterior samples
#' 
#' @param MCMCoutputs a list of r elements, where r is the number of EOFs specified, 
#'         and each element is a list of 6 elements containing posterior quantities of 
#'         parameters, log-likelihood, chisq statistics, and prior, which are 
#'         the outputs from MCMC algorithm
#' 
#' @param N an integer containing the number of ensemble memembers
#' 
#' @param n an integer containing the number of grid cells 
#' 
#' @param dir a file path specifying where to save the plots
#' 
#' @param burnin an integer specifying where the posterior samples have good mixing 
#' propertity
#' 
#' @param format a character specifying the format of plots to be saved, e.g., "pdf", "png"
#' 
#' @param width the width of the device (inches for .pdf plot and pixels for others)
#' 
#' @param height the height of the device (inches for .pdf plot and pixels for others)
#' 
#' @param weights_plot a logical value specifying whether the posterior plots of
#' weights are saved
#' 
#' @param beta_plot a logical value specifying whether the posterior plots of beta
#' are saved
#' 
#' @param beta_box a logical value specifying whether the box plots of posterior 
#' samples for beta are saved
#' 
#' @param chisq_plot a logical value specifying whether the chisquare test plot 
#' is saved
#' 
#' @return null
#' 
#' @description This function mainly saves relevant plots for posterior samples.
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export
#' 
#' @keywords figures
#' 
#' @seealso parDAbayes




plotMCMC <- function(MCMCoutputs, N, n, dir=getwd(), burnin=2000, format="pdf",
	width=8, height=6, weights_plot=TRUE, beta_plot=TRUE, beta_box=TRUE, chisq_plot=TRUE){

if (format=="pdf"){ # the size is in inches
  if (is.null(width) & is.null(height)){
    my_width <- 8
    my_height <- 6
  }else{
    my_width <- width
    my_height <- height
  }
}else{ # the size is in pixels
  if (is.null(width) & is.null(height)){
    my_width <- 480*2
    my_height <- 480*1.5
  }else{
    my_width <- width
    my_height <- height
  }
}

Res <- MCMCoutputs
r_num <- length(Res)
niter <- length(Res[[1]]$logsigma)
p <- dim(Res[[1]]$beta)[1]

loglik_conv <- matrix(NA, niter-burnin, r_num)
prior_conv <- matrix(NA, niter-burnin, r_num)
beta_conv <- array(NA, dim=c(p, niter-burnin, r_num))
for (r_ind in 1:r_num){
  loglik_conv[ ,r_ind] <- Res[[r_ind]]$loglik[(burnin+1):niter]
  prior_conv[ ,r_ind] <- Res[[r_ind]]$prior[(burnin+1):niter]
  beta_conv[ , ,r_ind] <- Res[[r_ind]]$beta[ ,(burnin+1):niter]
}
const <- max(loglik_conv)

lik_divByConst <- exp(loglik_conv - const)
marg_liks <- apply(lik_divByConst*prior_conv,2,mean)
weights <- marg_liks / sum(marg_liks)

######################################################################
# plot of the different r (number of EOFS) and their posterior weights
if (weights_plot){
	do.call(format, list(paste(dir, "/posterior_weights.", format, sep=""), 
	                     height=my_height, width=my_width))		
	plot(1:r_num, weights, xlab="Number of EOFs", ylab="Posterior weights",
	 pch=19, cex=0.6, cex.lab=0.8, cex.axis=0.8, xaxt='n')
	axis(1, at=seq(1, r_num+2, by=2), cex.axis=0.6, las=1)
	dev.off()
}


# calculate posterior mean of beta
mean_post_beta <- rep(NA, p)
post_beta <- matrix(NA, niter-burnin, p)
for (i in 1:p){
	post_beta[ ,i] <- beta_conv[i, , ]%*%weights
	mean_post_beta[i] <- mean(post_beta[ ,i])
}

xmax_density <- (max(beta_conv)+min(beta_conv))/2+1.4*(max(beta_conv) - min(beta_conv))
xmin_density <- (max(beta_conv)+min(beta_conv))/2-1.4*(max(beta_conv) - min(beta_conv))

# density plot of posterior distributions for beta 
if(beta_plot){

	# under different number of EOFs
	do.call(format, list(paste(dir, "/density_beta_individual.", format, sep=""), 
	                     height=my_height, width=my_width))		
	plot(density(beta_conv[1, ,1], bw=.01), type="l", xlab=expression(beta), lwd=1, col="red", 
		main=" ", xlim=c(xmin_density,xmax_density), 
		ylim=c(0, ceiling(max(density(beta_conv[1, ,1], bw=.01)$y))+2),
		cex=0.6, cex.axis=0.8, cex.lab=0.8)
	lines(density(beta_conv[2, ,1], bw=.01), type="l", lwd=1, col="green")
	for (r_ind in 2:r_num){
		lines(density(beta_conv[1, ,r_ind], bw=.01), type="l", lwd=1, col="red")
		lines(density(beta_conv[2, ,r_ind], bw=.01), type="l", lwd=1, col="green")
	}
	legend("topright", legend=c("anthropogenic", "natural"), 
		col=c("red","green"), lty=c(1,1), cex=0.8)
	dev.off()

	# BMA
	do.call(format, list(paste(dir, "/density_beta_BMA.", format, sep=""), 
	                     height=my_height, width=my_width))		
	plot(density(post_beta[ ,1], bw=.1), type="l", xlab=expression(beta), lwd=2, col="red", 
		main=" ", xlim=c(xmin_density,xmax_density), cex=0.6, cex.axis=0.8, cex.lab=0.8)
	lines(density(post_beta[ ,2], bw=.1), type="l", lwd=2, col="green")
	legend("topright", legend=c("anthropogenic", "natural"), 
		col=c("red","green"), lty=c(1,1), cex=0.8)
	dev.off()
}


# boxplot of posterior distribution
if(beta_box){
	beta_conv_ant <- beta_conv[1, , ]
	do.call(format, list(paste(dir, "/box_anthro_beta.", format, sep=""), 
	                     height=my_height, width=my_width))	
	boxplot(beta_conv_ant, xlab="Number of EOFs", ylab="anthropogenic forcing", 
		cex=0.6, cex.axis=0.8, cex.lab=0.8, pars=list(boxwex=0.5), 
		medcol="red", medlwd=1)
	dev.off()

	do.call(format, list(paste(dir, "/box_nat_beta.", format, sep=""), 
	                     height=my_height, width=my_width))	
	beta_conv_nat <- beta_conv[2, , ]
	boxplot(beta_conv_nat, xlab="Number of EOFs", ylab="natural forcing", 
		cex=0.6, cex.axis=0.8, cex.lab=0.8, pars=list(boxwex=0.4), 
		medcol="red", medlwd=1)
	dev.off()
}


# residual consistency test
chisq_conv <- matrix(NA, niter-burnin, r_num)
for (r_ind in 1:r_num){
	chisq_conv[ ,r_ind] <- Res[[r_ind]]$chisq[niter-burnin]
}

pvals <- pchisq(c(chisq_conv), N*n)
pvals_2sided <- 1 - abs(pvals-0.5)*2
weights_all <- rep(weights, each=niter-burnin)

if(chisq_plot){
  normalized_weights <- exp(log(weights_all)-log(sum(weights_all)))
  pvalueData <- density(pvals_2sided, bw=.01, weights=normalized_weights)
  flag <- pvalueData$x > 0.990 
  pvalueData$x[flag==TRUE] <- 1
  pvalueData$y[flag==TRUE] <- 0
  pvalmin <- min(pvalueData$x)
  pvalmax <- max(pvalueData$x)
  do.call(format, list(paste(dir, "/chisq_consistency_test.", format, sep=""), 
	                     height=my_height, width=my_width))
	plot(pvalueData$x, pvalueData$y, main=" ", xlab="p-value", ylab="density", 
		type="l", xlim=c(pvalmin, pvalmax), cex=0.6, cex.axis=0.8, cex.lab=0.8)
	dev.off()
}



}
