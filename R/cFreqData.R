#' @title Creates a cFreqData class object
#' @description Creates a \code{cFreqData} class object from user supplied time frequency data.
#' @param x A matrix storing the frequency of each cultural variant (columns) at each sampling phase (rows).
#' @param timestamp A vector containing the timestamp (i.e. last timestep) of each sampling phase.
#' @param duration A vector containing the durations of the sampling phase (except the first phase) calculated in number of timesteps. If a single value is supplied the durations will be assumed to be the same for all phases. 
#' @return An object of class cFreqData.
#' @examples
#' dt=as.matrix(data.frame(var1=c(0,10,15),var2=c(20,18,0),var3=c(0,0,4),var4=c(32,20,18)))
#' x=cFreqData(dt,timestamp = c(1,6,11))
#' @import stats
#' @import utils 
#' @export

cFreqData <- function(x,timestamp,duration=1){
  
  # Check matrix and length of starttimes
  if (nrow(x) != length(timestamp))
  {
    stop("The number of rows in the frequency matrix does not correspond to the number of timestamps")
  }
  
  # Check matrix and length of starttimes
  if (length(duration)!=c(nrow(x)-1) & length(duration)!=1)
      {
  stop("The vector duration should have length of 1 or the number of sampling phases minus 1")  
      }

  if (length(duration)==1) 
  {
	duration=rep(duration,nrow(x)-1)
  }

  # Convert timestep for convenience
  timestamp=timestamp-timestamp[1]+1
 
  # Check if durations make sense
  if (!all(diff(timestamp)>=duration))
  {
  stop("Supplied phase duration are too long")  	
  }


  # Reorder matrix:

    if (any(apply(x,2,sum)==0)){x=x[,-which(apply(x,2,sum)==0)]}
    x = x[,order(x[1,],decreasing=TRUE)]	

    tmp=which(x[1,]==0)
    soFar=min(tmp)-1

    for (i in 2:c(nrow(x)))
        {
            if (any(x[i,tmp]>0))
                {
                    count=length(which(x[i,tmp]>0))
                    x=x[,c(1:soFar,order(x[i,tmp],decreasing=TRUE)+soFar)]
                    soFar= soFar+ count	
                    tmp=which(apply(x[1:i,],2,sum)==0)
                }
        }

  # Compute Target frequencies
  freq=prop.table(x,1) 
  freq=freq[-1,which(freq[1,]>0)]
  target.freq=as.numeric(freq)

  
  res=list(cfreq=x,k=apply(x,1,function(x){sum(x>0)}),n=apply(x,1,sum),n.phases=nrow(x),tstamp=timestamp,duration=duration,target.freq=target.freq)
  class(res) <- c("cFreqData",class(res))
  return(res)
}
