#' @title Plot posterior predicted checks
#'
#' @description Plots pairwise joint and marginal distributions of ABC posterior samples.  
#' @param x a matrix containing the posterior samples.
#' @param pnames a vector containing parameter names
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export  


plotPost<-function(x,pnames)
{
	n.param = ncol(x)
	par(mfcol=c(n.param,n.param))
	par(mar=c(7,5,1,1))
	nsim = nrow(x)
	if (nsim<1000){nsim=1000}

	quadSize = n.param*n.param
	position = matrix(1:quadSize,n.param,n.param)
	marginal = diag(position)
	joint = which(lower.tri(position))
	outer = which(upper.tri(position))
	combinations = expand.grid(1:n.param,1:n.param)

	for (i in 1:quadSize)
	{
	  if (i%in%marginal)
	  {
	    k = which(marginal==i)
	    tmp=density(x[,k])
	    plot(0,0,xlim=range(tmp$x),ylim=range(tmp$y),type="n",axes=F,xlab="",ylab="")
	    polygon(c(tmp$x,rev(tmp$x)),c(tmp$y,rep(0,length(tmp$y))),col="deepskyblue3",border=NA)
	    abline(v=median(x[,k]),lty=2,lwd=1.5)
	    axis(1)
	    axis(2)
	    mtext(pnames[k],1,line=3)
	  }
	  if (i%in%joint)
	  {
	    v = combinations[i,1]
	    u = combinations[i,2]
	    
	    smoothScatter(x[,u],x[,v],xlab="",ylab="")
	    axis(1)
	    axis(2)
	    mtext(pnames[u],1,line=3)
	    mtext(pnames[v],2,line=3,las=2)
	  }
	  if (i%in%outer)
	  {
	    plot(0,0,type='n',axes=F,xlab='',ylab='')
	  }
	}
	par(mfrow=c(1,1))
	
	  
	  
	}
	
	
