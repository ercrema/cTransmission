#' @title Plot posterior predicted checks
#'
#' @description Compares observed variants frequency to the expected variants frequency of the inferred transmission model. 
#' @param x A \code{pcheck} class object obtained using \code{link{predCheck}}.
#' @param var.names Variant names supplied with the same order as the observed data.
#' @param phase.names Names of the sampling phases
#' @param index Integer value indicating which sampling phase to display counting from the second observed sampling phase.
#' @param ppi Prediction percentile interval. Default is 0.95.
#' @seealso \code{\link{predCheck}}
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @import LaplacesDemon
#' @export  


plot.pcheck <- function(x,var.names=NULL,phase.names=NULL,index=1,ppi=0.95,...)
{
	xrange = range(min(c(x$observed,unlist(x$predicted))),max(c(x$observed,unlist(x$predicted))))	
	
	if (is.null(var.names))
	{
		var.names=x$var.names
		if (is.null(var.names))
		{
			warning("Variable names are not supplied, using column number instead")
			var.names = 1:ncol(x$observed)
		}
	}

	if (is.null(phase.names))
	{
		phase.names=x$phase.names[-1]
		if (is.null(phase.names))
		{
			warning("Phase names are not supplied, using row number instead")
			phase.names = 1 + 1:nrow(x$observed)
		}
	}

# plot window	
	plot(0,0,type="n",xlim=xrange+c(-0.1,0.1),ylim=c(0.5,ncol(x$observed)+0.5),axes=F,xlab="Relative Frequencies",ylab="Variants",main=paste0("Phase ",phase.names[index]),...)

	axis(side=1)
	axis(side=2,at=1:ncol(x$observed),labels=var.names,las=2)
	box()
# compute credible interval

	pred=x$predicted[[index]]
	colnames(pred)=var.names
	hinterval = t(apply(pred,2,quantile,prob=c((1-ppi)/2,ppi+(1-ppi)/2)))
# 	hinterval = LaplacesDemon::p.interval(pred,prob=hpdi)
	
# plot credible intervals and observed frequencies

	for (i in 1:ncol(pred))
	{
	   rng = hinterval[i,] 
	   arrows(x0=rng[1],x1=rng[2],y0=i,y1=i,code=3,length=0.05,angle=90)
	   colpts = "black"
	   if (x$observed[index,i]>rng[2]|x$observed[index,i]<rng[1])
	   {
		   colpts="red"
	   }
	   points(x=x$observed[index,i],y=i,pch=20,col=colpts,cex=1.5)
	}
}

