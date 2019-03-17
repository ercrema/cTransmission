#' @title Frequency Biased Transmission Model
#' @description Single parameter frequency biased transmission model
#' @param x Relative or absolute frequency of cultural variants.
#' @param k Number of initial cultural variants with the last slot for mutation events (see details below). 
#' @param mu Innovation rate.
#' @param b Transmission bias prameter. When b>0 the transmission has conformist bias, when b<0 the transmission has anticonformist bias, and when b=0 the transmission is unbiased.
#' @details. The function estimates the probability of  each variant being selected in the tranmission process. The last slot of the user supplied vector \code{x} is designated to the outcome of innovation events.
#' @export


frequencyBias<-function(x,k,mu,b)
{
    h = rep(0,k)
    transProb = c(0,(k))
    x = as.vector(x / sum(x))
    h = x^(1 + b)
    
    #Normalisation
    h = h / sum(h)  
	  transProb[1:k-1] = h[1:k-1] * (1 - mu)
	  transProb[k] = h[k] + (sum(h[1:k-1])) * mu
	  
	  return(transProb)	
}


#' @title Conformist-Biased Transmission Model
#' @param x Relative or absolute frequency of cultural variants.
#' @param k Number of initial cultural variants with the last slot for mutation events (see details below). 
#' @param mu Innovation rate.
#' @param c Probability of an inividual engaging into a conformist biased transmission (i.e. the selection of the most common cultural variant). 
#' @details. The function estimates the probability of  each variant being selected in the tranmission process. The last slot of the user supplied vector \code{x} is designated to the outcome of innovation events.

conformistBias<-function(x,k,mu,c)
{
    h = rep(0,k)
    transProb = c(0,(k))
    x = as.vector(x / sum(x))

    h = (1-c)*x
    h[1] = h[1]+c
    
    
    #Normalisation
    transProb[1:k-1] = h[1:k-1] * (1 - mu)
    transProb[k] = h[k] + (sum(h[1:k-1])) * mu
	  
    return(transProb)	
}




