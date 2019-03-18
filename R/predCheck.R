#' @title Posterior predictive check for cultural transmission model.
#' @description Creates a \code{pcheck} class object containing expected variant frequencies based on posterior parameter values.
#' @param x A cFreqData class object. 
#' @param sim.model A simulation model created using the \code{genSim} function.
#' @param posterior Posterior samples. Can be extracted from the output of \code{EasyABC} functions such as \code{ABC_rejection} and \code{ABC_sequential}. 
#' @param nsim Number of simulations for posterior predictive check. 
#' @return An object of class pcheck.
#' @seealso \code{\link{genSim}}; \code{\link{cFreqData}}
#' @import stats
#' @import utils 
#' @export

predCheck <- function(x,sim.model,posterior,nsim=NULL)
{
	if (nsim > nrow(posterior)|is.null(nsim))
	{
		nsim = nrow(posterior) 
		warning(paste0("Number of simulations undefined or higher than the posterior sample size. nsim has been set to ", nsim))
	}
	
	# create placeholder
	predicted = vector("list", length=x$n.phases-1)
	for (i in 1:length(predicted)){predicted[[i]]=matrix(NA,nrow=nsim,ncol=x$k[1])}

	# make simulations
	posterior = posterior[sample(1:nrow(posterior)),]
	for (i in 1:nsim)
	{
		tmp = sim.model(as.list(posterior[i,]))
	        tmp = matrix(tmp,nrow=x$n.phases-1,ncol=x$k[1])
		for (p in 1:nrow(tmp))
		{
			predicted[[p]][i,]=tmp[p,]
		}
	}
	
	# store output
	observed = matrix(x$target.freq,nrow=x$n.phases-1,ncol=x$k[1])
	var.names = colnames(x$cfreq)[1:x$k[1]]
	phase.names = rownames(x$cfreq)	
	result=list(observed=observed,predicted=predicted,var.names=var.names,phase.names=phase.names)
	class(result) = c("pcheck",class(result))
	return(result)
}


