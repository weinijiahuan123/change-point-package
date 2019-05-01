#' Get Change Points from Peak Ranges
#'
#' Transform the peak ranges of change points to exact change points.
#'
#' Find the exact change points with peak ranges based on log likelihood comparison.
#' @param y The original data to find change points. Must be one dimensional data.
#' @param n_peak_range The number of peak ranges of change points
#' @param peak_range The location of ranges of change points
#' @param get_loglik The method to get
#'
#' @return change_point
#' @export
#' @examples
#' N <- 1000
#' N1 <- floor(0.1*N)
#' N2 <- floor(0.3*N)
#' a1 <- c(0.8, -0.3); c1 <- 0
#' a2 <- c(-0.5, 0.1); c2 <- 0
#' a3 <- c(0.5, -0.5); c3 <- 0
#' y <- rep(0,N)
#' L<-2
#' y[1:L] <- rnorm(L)
#' for (n in (L+1):N){
#'   if (n <= N1) {
#'     y[n] <- y[(n-1):(n-L)] %*% a1 + c1 + rnorm(1)
#'   } else if (n <= (N1+N2)) {
#'     y[n] <- y[(n-1):(n-L)] %*% a2 + c2 + rnorm(1)
#'   }
#'   else {
#'     y[n] <- y[(n-1):(n-L)] %*% a3 + c3 + rnorm(1)
#'   }
#' }
#' RangeToPoint(y,n_peak_range=2,peak_range=list(seq(70,105),seq(395,420)))
RangeToPoint<-function(y,n_peak_range,peak_range,get_loglik=GetLogLik) {
  len <- length(y)
  # Initialize change points
  change_point <- numeric(n_peak_range)
  for (i in 1:n_peak_range) {
    change_point[i] <- ceiling(mean(peak_range[[i]]))
  }
  # add 1 as the first change point and N as the last change point for computation convenience
  change_point <- c(1, change_point, N)

  for (i in 1:n_peak_range) {
    best_log_lik <- Inf
    #test
    #i=1
    #test
    #cat("i=",i,"\n_peak_range")
    # Get the outer range to find the exact change point.

#
#     if (i == n_peak_range) {
#       outer_range<-y[change_point[i]:len]
#     } else {
#       outer_range<-y[change_point[i]:(peak_range[[i+1]][1]-1)]
#     }

    # Fix following ranges and find the exact change point of range i with smallest quadratic loss.
    for (e in 1:length(peak_range[[i]])) {
      #test
      #e=1
      #test
      left_part <- get_loglik(y, change_point[i]+1, peak_range[[i]][e])
      right_part <- get_loglik(y, peak_range[[i]][e]+1, change_point[i+2])
      log_lik <- left_part+right_part

      if (log_lik < best_log_lik) {
        best_log_lik <- log_lik
        change_point[i+1] <- peak_range[[i]][e]
      }
      #cat("e",e,"\n")
      #cat("loss",loss,"\n")
      #cat("best_log_lik",best_log_lik,"\n")
      #cat("change_point",change_point,"\n")
    }
  }
  # Delete the first and last change points
  change_point <- change_point[-c(1,length(change_point))]
  return(change_point)
}


#' Get Log Likelihood
#'
#' For a series of one dimensional data, get the log likelihood.
#' @param y The data to calculate log likelihood. The data must be one dimesional.
#' @param left The left end of the data that users want to use.
#' @param right The right end of the data that users want to use.
#'
#' @return log_lik
#' @export
#'
#' @examples
GetLogLik <- function(y, left, right) {
  EstimateAr <- function(x,T1,T2,L){
    if (T1 > (T2 - L)) {
      warning("Error in estimate_ar")
    }
    if (T1 <= L) {
      T1 <- L + 1
    }
    Y <- matrix(0, nrow=L+1, ncol=T2-T1+1)
    Y[1,] <- 1
    for (k in 1:L) {
      Y[k+1, ] <- x[(T1-k):(T2-k)]
    }
    A <- Y%*%t(Y)
    B <- Y%*%x[T1:T2]
    C <- solve(A)%*%B
    e <- x[T1:T2]-t(Y)%*%C
    sigma2 <- sum(e^2)/(T2-T1+1)
    est_coef <- list(C=C,sigma2=sigma2)
    return(est_coef)
  }

  log_lik <- (right-left)*log(EstimateAr(y,left+1,right,2)$sigma2)
}
