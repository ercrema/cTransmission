#' @import stats
#' @keywords internal

instances<-function(x,variants)    
{
    x=c(x,variants)
    res<-table(x)-1
    return(res)    
}


