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

	for (i in 1:n.param)
	{
		for (j in 1:n.param)
		{
			# Marginal Distribution
			if (i==j)
			{
				tmp=density(x[,i])
				plot(0,0,xlim=range(tmp$x),ylim=range(tmp$y),type="n",axes=F,xlab="",ylab="")
				polygon(c(tmp$x,rev(tmp$x)),c(tmp$y,rep(0,length(tmp$y))),col="deepskyblue3",border=NA)
				abline(v=median(x[,i]),lty=2,lwd=1.5)

				axis(1)
				axis(2)
				mtext(pnames[i],1,line=3)
			}

			# Joint Distribution
			if (i<j)
			{
				smoothScatter(x[,i],x[,j],xlab="",ylab="")
				axis(1)
				axis(2)
				mtext(pnames[i],1,line=3)
				mtext(pnames[j],2,line=3,las=2)
			}

			if (j>i)
			{
			        plot(0,0,type='n',axes=F,xlab='',ylab='')
			}

		}

	}
	par(mfrow=c(1,1))
}

