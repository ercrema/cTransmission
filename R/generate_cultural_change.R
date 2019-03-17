#' @title Computes exected frequency of cultural variants. 
#' @description This functions computes expected frequencies of extant cultural variants given a user defined transmission model and observed data.
#' @param x A cFreqData class object.
#' @param iniPop Raw population frequency of culturial variants
#' @param rr Number of individuals to be removed and added at each tranmission even. Output of \code{generate_removalReplacement}.
#' @param mu Innovation rate.
#' @param param A list containing parameters values to be forwarded to the transmission model.
#' @param model A transmission model (see details below)
#' @details Transmission model can be either the one provided by \code{cTransmission} or custom generated. Custom generated model should be written as a function which outputs the probability of selecting each variant. The function would require at least the arguments x (the raw frequencies of each variant), k (the number of variants), and mu (the mutation rate), plus additional parameters. 
#' @import stats
#' @import utils 
#' @keywords internal



generate_cultural_change<-function(x, iniPop, rr, mu, params=list(b=0),model=frequencyBias)
{
  # Extract variables from observed data
    nSamples = x$n.phases #length(samples) #Number of samples
    nType = x$k[1]+1 #length(samples[[1]])+1 #Number of variant types in the first phase
    timePoints = x$tstamp
    durations = x$duration
    u = rr$u
    v = rr$v
    sampleSize = x$n
    	


  # Create matrix storing all variants
    timesteps = timePoints[1]:timePoints[length(timePoints)]  	
    nts = length(timesteps)
    resmat= matrix(0,nrow=nts,ncol=nType)
    resmat[1,]=iniPop

  # Matrix storing the sampled variants
    sampmat = matrix(NA,nrow=nSamples-1,ncol=nType) 



    for (t in 1:c(nts-1))
    {

	    relFreq = resmat[t,] / sum(resmat[t,])
	    #Calculation of transmission probabilities
# 	    transProb = generate_transmissionProb(nType,resmat[t,],b,mu)

	    #flexible approach
	    arguments = c(list(x=resmat[t,],mu=mu,k=nType),params)
	    transProb = do.call(model,arguments)
        

	    #Removal of u[t] variants. First a check:
	    if (u[t]>sum(resmat[t,]))
	    {
		    stop("Number of variants to be removed is larger than population size")
	    }	      

	    if (u[t]<=sum(resmat[t,]))
	    {

		    h = sample(1:nType,u[t], replace = TRUE, prob = relFreq) 
		    h = instances(h,1:nType)
		    resmat[t+1,]= resmat[t,] - h

		    while (any(resmat[t+1,]<0))
		    {
			    resmat[t+1,]=NA
			    h = sample(1:nType,u[t], replace = TRUE, prob = relFreq) 
			    h = instances(h,1:nType)
			    resmat[t+1,]= resmat[t,] - h
		    }
	    }

	    #Addition of v[t] variants according to probability transProb
	    h = sample(1:nType, v[t], replace = TRUE, prob = transProb)
	    h = instances(h,1:nType)
	    resmat[t+1,]= resmat[t+1,] +  h
    }

    #Sampling based on duration
    for (i in 1:c(nSamples-1))
    {
	st = timePoints[i+1] - durations[i] + 1
	en = timePoints[i+1]    
	tmp=resmat[st:en,]
	tmp=apply(tmp,2,sum)
	h=sample(1:nType,size=sampleSize[i+1],replace=TRUE,prob=tmp/sum(tmp))
	sampmat[i,]=instances(h,1:nType)
    }

   #compute relative frequencies 
   theosamples = prop.table(sampmat,1)
   # remove the mutation column
   theosamples = as.numeric(theosamples[,-nType])   
   return(theosamples)
  }  
	

	
