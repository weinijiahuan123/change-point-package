RangeToPoint<-function(y,n_peak_range,peak_range,get_loglik=EstimateAr) {
  # Transform the peak ranges of change points to exact change points.
  #
  # Args:
  #   y: The original data to find change points. Must be one dimensional data
  #   n_peak_range: The number of peak ranges of change points
  #   peak_range: The location of ranges of change points
  #
  # Returns:
  #   change_point: The location of exact change points
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
      left_part <- (peak_range[[i]][e]-change_point[i])*log(get_loglik(y,change_point[i]+1,peak_range[[i]][e],2)$sigma2)
      right_part <- (change_point[i+2]-peak_range[[i]][e])*log(get_loglik(y,peak_range[[i]][e]+1,change_point[i+2],2)$sigma2)
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
a=seq(70,105)
b=seq(395,420)
peak_range<-list(a,b)
RangeToPoint(y,n_peak_range,peak_range)
res<- MultiWindow(y,window_list=c(100,50,20,10,5),point_max=5,seg_min=2,tolerance=1)
n_peak_range=res$n_peak_range
peak_range=res$peak_range

EstimateAr=function(x,T1,T2,L){
  if (T1>(T2-L)) {
    warning("Error in estimate_ar")
  }
  if (T1<=L) {
    T1=L+1
  }
  Y=matrix(0,nrow=L+1,ncol=T2-T1+1)
  Y[1,]=1
  for (k in 1:L) {
    Y[k+1,]=x[(T1-k):(T2-k)]
  }
  A=Y%*%t(Y)
  B=Y%*%x[T1:T2]
  C=solve(A)%*%B
  e=x[T1:T2]-t(Y)%*%C
  sigma2=sum(e^2)/(T2-T1+1)
  est_coef=list(C=C,sigma2=sigma2)
  return(est_coef)
}
