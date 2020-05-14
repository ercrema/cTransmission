#' @title Approximate Bayesian Computation inference via Rejection Method
#' @description Estimates model parameter via Approximate Bayesian Computation.
#' @param x A cFreqData class object with the observed variant frequencies.
#' @param sim.model A transmission model generated using the \code{genSim()} function.
#' @param prior A list of model priors. Each element of the list corresponds to a model parameter. The list element must be a vector whose first argument determines the type of prior distribution: possible values are "unif" for a uniform distribution on a segment and "normal" for a normal distribution. The following arguments of the list elements contain the characteritiscs of the prior distribution chosen: for "unif", two numbers must be given: the minimum and maximum values of the uniform distribution; for "normal", two numbers must be given: the mean and standard deviation of the normal distribution. Syntax adapted from the \code{EasyABC} package.
#' @param nsim Number of simulations
#' @param tol Required proportion of accepted simulations.
#' @param raw Determines whether raw simulation output and parameters should be returned. Default is FALSE.
#' @param ncore Number of cores for parallel processing. Default is 1.
#' @details 
#' @import utils
#' @import abc
#' @import stats
#' @import foreach
#' @import parallel
#' @import doParallel
#' @export

abcRej<-function(x,sim.model,prior,nsim,tol,ncore=1,raw=FALSE)
{

	#Construct prior
	nparam=length(prior)
	param=matrix(NA,nrow=nsim,ncol=nparam)


	for (i in 1:nparam)
	{
		model=prior[[i]][1]
		if (model=='unif')
		{
			param[,i]=runif(nsim,as.numeric(prior[[i]][2]),as.numeric(prior[[i]][3]))
		}
		if (model=='norm')
		{
			param[,i]=rnorm(nsim,as.numeric(prior[[i]][2]),as.numeric(prior[[i]][3]))
		}
	}

	# Execute simulation

	cl <- makeCluster(ncore)
	registerDoParallel(cl)

	simres<-foreach (i=1:nsim,.combine=rbind,.packages='cTransmission') %dopar% {
		p = as.list(param[i,])
		sim.model(p)
	}


	#ABC rejection inference
	result=abc(x$target.freq,param,simres,tol=tol,method='rejection')
	if (!raw){return(result)}
	if (raw){return(list(abcRes=result,simres=simres,param=param))}
}











