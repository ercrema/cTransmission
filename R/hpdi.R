#' @import stats
#' @keywords internal


hpdi <- function (obj,prob = 0.95) 
{
 
    vals <- apply(obj, 2, sort)
    if (!is.matrix(vals)) 
      stop("obj must have nsamp > 1.")
    
    nsamp <- nrow(vals)
    npar <- ncol(vals)
    gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
    init <- 1:(nsamp - gap)
    inds <- apply(vals[init + gap, , drop = FALSE] - vals[init,, drop = FALSE], 2, which.min)
    ans <- cbind(vals[cbind(inds, 1:npar)], vals[cbind(inds + gap, 1:npar)])
    
    rownames(ans) <- colnames(obj)
      
    vals <- apply(obj, 2, sort)
      if (!is.matrix(vals)) 
        stop("obj must have nsamp > 1.")
      for (m in 1:ncol(vals)) {
        kde <- density(vals[, m])
        dens <- approx(kde$x, kde$y, vals[, m])$y
        dens.ind <- dens >= as.vector(quantile(dens, probs = 1 - prob)) * 1
        ints <- ""
        count <- 1
        for (i in 1:nrow(vals)) {
          if ((i == 1) & (dens.ind[i] == 1)) {
            ints <- paste("(", round(vals[i, m], 3), 
                          ",", sep = "")
            if (count > ncol(ans)) 
              ans <- cbind(ans, NA)
            ans[m, count] <- vals[i, m]
            count <- count + 1
          }
          if (i > 1) {
            if ((dens.ind[i] == 0) & (dens.ind[i - 1] ==1)) {
              ints <- paste(ints, round(vals[i - 1, m], 3), ")", sep = "")
              if (count > ncol(ans)) 
                ans <- cbind(ans, NA)
              ans[m, count] <- vals[i - 1, m]
              count <- count + 1
            }
            if ((dens.ind[i] == 1) & (dens.ind[i - 1] == 0)) {
              ints <- paste(ints, " (", round(vals[i, m], 3), ",", sep = "")
              if (count > ncol(ans)) 
                ans <- cbind(ans, NA)
              ans[m, count] <- vals[i, m]
              count <- count + 2
            }
          }
        }
        if ((dens.ind[i] == 1) & (dens.ind[i - 1] == 1)) {
          ints <- paste(ints, round(vals[i, m], 3), ")",  sep = "")
          if (count > ncol(ans)) 
          ans <- cbind(ans, NA)
          ans[m, count] <- vals[i, m]
          count <- count + 1
        }
      }
    
  return(ans)
}
