#' @title Generates a simulation model for Approximate Bayesian Computation inference.
#' @description Generates an R function that can be used for generative inference.
#' @param theta A vector of parameter names required for the transmission model and to be estimated using ABC (see details below).
#' @param x A cFreqData class object containing the observed frequency of cultural variants.
#' @param model A transmission model (e.g. \code{frequencyBias}). See vignette for more details.
#' @param alpha Shape parameter for the Dirichlet distribution required for estimating population relative frequencies.
#' @param sMean Mean sampling fraction (if not supplied as argument in \code{theta}).
#' @param sVariance Variance (uncertainty) of sampling fraction
#' @param rMean mean fraction of population being replaced at each timestep (if not supplied as argument in \code{theta}).
#' @param rVariance Variance (uncertainty) of fraction of population being replaced at each timestep
#' @param popSize Vector of population sizes for each sampling phase. Required if sampling fraction is not specified elsewhere.
#' @details The argument \code{theta} should contain the name of model parameters to be inferred using Approximate Bayesian Computation. The names "r", "s", and "mu" refer specifically to replacement rate, sampling fraction, and innovation rate. All other names will be supplied to the transmission model. Notice that when executing the transmission model or using this in ABC inference the order of with which the values are supplied should match the same order provided for \code{theta}. 
#' @examples
#' \donotrun{
#' dt=as.matrix(data.frame(var2=c(0,23,50),var2=c(200,180,180),var3=c(0,10,40),var4=c(320,290,270)))
#' x=cFreqData(dt,timestamp = c(1,6,11), duration=c(2,3))
#' sim = genSim(theta=c("s","mu","b"),x=x,model=frequencyBias,alpha=1,rMean=0.2,rVariance=0.01)
#' sim(list(0.2,0.01,0))
#' prior=list(prior_s=c("unif",0.3,0.5),prior_mu=c("unif",0.001,0.01),prior_b=c("normal",0,0.1))
#' res = abcRej(x=x,sim.model=sim,prior=prior,nsim=1000,tol=0.1,ncore=2)
#' ## Using the EasyABC Package
#' library(EasyABC)
#' res = ABC_rejection(sim,prior=list(prior_s=c("unif",0.3,0.5),prior_mu=c("unif",0.001,0.01),prior_b=c("normal",0,0.1)),tol=0.01,nb_simul=1000) 
#' }
#' @export


genSim<-function(theta, x, model, mu=NULL, sMean=NULL, sVariance=NULL, alpha, rMean = NULL, rVariance = NULL, popSize = NULL)
{
 # Extract key variables
  nSamples = x$n.phases 
  nType = x$k[1]+1 
  iniSample = c(x$cfreq[1,],0) 
  sampleSize=x$n 

# Extract parameters from theta
  r.index = s.index = mu.index = NULL
  if (any(theta=="r"))
  {
	  r.index = which(theta=="r")
	  if (!is.null(rMean)|!is.null(rVariance))
	  {
		stop("Replacement rate r cannot be both a fixed value and a prior")
	  }
  }  
  if (any(theta=="s"))
  {
	  s.index = which(theta=="s")
	  if (!is.null(sMean)|!is.null(sVariance))
	  {
		stop("Sampling fraction s cannot be both a fixed value and a prior")
	  }
  }  
  if (any(theta=="mu"))
  {

	  mu.index = which(theta=="mu")
	  if (!is.null(mu))
	  {
		  stop("Innovation rate mu cannot be both a fixed value and a prior")
	  }

  }
  if (any(!theta%in%c("r","s","mu"))){transmission.index = which(!theta%in%c("r","s","mu"))}


# function to be forwarded to EasyABC  
  tSamples<-function(params){
    require(cTransmission)
# Extract r and s if present as prior:
  if (!is.null(r.index)) {r=params[[r.index]]}	  
  if (!is.null(s.index)){sMean=params[[s.index]];sVariance=0}

# Infer population sizes:
  if(is.null(popSize)&(is.null(sMean)|is.null(sVariance)))
  {
	  stop("Either popSize or sMean & sVariance should be provided")
  }

  if(is.null(popSize)){
	  popSize = generate_popSize(x, sMean, sVariance)
  }

  if (length(popSize) != x$tstamp[x$n.phases]-x$tstamp[1]+1){
	  print("error : popSize does not have the right size")
  }

# Infer variants distribution at t1  

  iniPop=round(generate_initialPop(c(x$cfreq[1,which(x$cfreq[1,]>0)],0),nsim=1,alpha = 1)*popSize[1])
  popSize[1] = sum(iniPop) #update popSize 1 in case discrepancies due to rounding
    
# Replacement Rate
  if (is.null(r.index))
  {
	  if (is.null(rMean)&is.null(rVariance))
	  {
		  stop("r or rMean and rVariance should be provided")
	  }
	  r = rnorm(1,rMean,rVariance)
	  while(r <= 0 | r >1){ #NOTE: ensure r is also not larger than 1?
		  r = rnorm(1,rMean,rVariance)
	  }
  }
  rr = generate_removalReplacement(x,N=popSize,r=r)    
 
  # Cultural Change
  if (!is.null(mu.index)){mu=params[[mu.index]]}
  if (is.null(mu.index)&is.null(mu))
  {
	  stop("Innovation rate mu should be provided")
  }

  #transmission parameters
  transmission.param=params[transmission.index]
  theoSamples = generate_cultural_change(x=x,rr=rr,iniPop=iniPop,params=transmission.param,mu=mu,model=model)
  return(theoSamples)
  }
  return(tSamples)
}








culturalChange_ABC_model.bkp<-function(theta, x, model, mu=NULL, sMean=NULL, sVariance=NULL, alpha, rMean = NULL, rVariance = NULL, popSize = NULL)
{
 # Extract key variables
  nSamples = x$n.phases 
  nType = x$k[1]+1 
  iniSample = c(x$cfreq[1,],0) 
  sampleSize=x$n 

# Extract parameters from theta
  r.index = s.index = mu.index = NULL
  if (any(names(theta)=="r")){r.index = which(names(theta)=="r")}  
  if (any(names(theta)=="s")){r.index = which(names(theta)=="s")}  
  if (any(names(theta)=="mu")){mu.index = which(names(theta)=="mu")}
  if (any(!names(theta)%in%c("r","s","mu"))){transmission.index = which(!names(theta)%in%c("r","s","mu"))}


# function to be forwarded to EasyABC  
  tSamples<-function(theta){


# Extract r and s if present as prior:
  if (!is.null(r.index)) {r=theta[[r.index]]}	  
  if (!is.null(s.index)){sMean=theta[[s.index]];sVariance=0}

# Infer population sizes:
  if(is.null(popSize)&(is.null(sMean)|is.null(sVariance)))
  {
	  stop("Either popSize or sMean & sVariance should be provided")
  }

  if(is.null(popSize)){
	  popSize = generate_popSize(x, sMean, sVariance)
  }

  if (length(popSize) != x$tstamp[x$n.phases]-x$tstamp[1]+1){
	  print("error : popSize does not have the right size")
  }

# Infer variants distribution at t1  

  iniPop=round(generate_initialPop(c(x$cfreq[1,which(x$cfreq[1,]>0)],0),nsim=1,alpha = 1)*popSize[1])
  popSize[1] = sum(iniPop) #update popSize 1 in case discrepancies due to rounding
    
# Replacement Rate
  if (is.null(r.index))
  {
	  if (is.null(rMean)&is.null(rVariance))
	  {
		  stop("r or rMean and rVariance should be provided")
	  }
	  r = rnorm(1,rMean,rVariance)
	  while(r <= 0 | r >1){ #NOTE: ensure r is also not larger than 1?
		  r = rnorm(1,rMean,rVariance)
	  }
  }
  rr = generate_removalReplacement(x,N=popSize,r=r)    
 
  # Cultural Change
  if (!is.null(mu.index)){mu=theta[[mu.index]]}
  if (is.null(mu.index)&is.null(mu))
  {
	  stop("Innovation rate mu should be provided")
  }

  transmission.param=theta[transmission.index]

  theoSamples = generate_cultural_change(x=x,rr=rr,iniPop=iniPop,params=transmission.param,mu=mu,model=model)
  return(theoSamples)
  }
  return(tSamples)
}







