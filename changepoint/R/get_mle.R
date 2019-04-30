#' Estimate Coefficients
#'
#' Transform N dependent data into approximated independent data (N/window_size) x (L+1).
#' Computes the estimated coefficients of each window of original data.
#'
#' @references
#' J. Ding, Y. Xiang, L. Shen, and V. Tarokh, \emph{Multiple Change Point Analysis:
#' Fast Implementation and Strong Consistency}. IEEE Transactions on Signal
#' Processing, vol. 65, no. 17, pp. 4495-4510, 2017.
#' @param y The original data to find change points.
#' @param window_size The number of observations each window contains.
#' @param L Lag order of the dataset. L>=1
#' @param method The method used to estimate coefficients.
#'           Must be one of c("yule-walker", "ols", "mle", "yw").
#'           Default is ols (ordinal least squares).
#'
#' @return x: The transformed data, which are the estimated coefficients of original data.
#' @export
#' @examples
#' N = 1000
#' N1 = floor(0.1*N)
#' N2 = floor(0.3*N)
#' a1 = c(0.8, -0.3); c1 = 0
#' a2 = c(-0.5, 0.1); c2 = 0
#' a3 = c(0.5, -0.5); c3 = 0
#' x = rep(0,N)
#' L=2
#' x[1:L] = rnorm(L)
#' for (n in (L+1):N){
#'   if (n <= N1) {
#'     x[n] = x[(n-1):(n-L)] %*% a1 + c1 + rnorm(1)
#'   } else if (n <= (N1+N2)) {
#'     x[n] = x[(n-1):(n-L)] %*% a2 + c2 + rnorm(1)
#'   }
#'   else {
#'     x[n] = x[(n-1):(n-L)] %*% a3 + c3 + rnorm(1)
#'   }
#' }
#' GetMle(x,window_size=100,L=2,method="mle")
#' GetMle(x,window_size=100,L=2,method="ols")
#' GetMle(x,window_size=100,L=2,method="yw")
GetMle <- function(y, window_size) {
  N <- length(y)
  n_window <- ceiling(N/window_size)
  L <- 2
  x <- matrix(0, nrow = n_window, ncol = L+1)
  for (n in 1:n_window) {
    #test
    #test
    #get estimated coefficients including constant
    est <- ar(y[(1 + (n - 1) * window_size):min(n * window_size, N)], aic = FALSE, order.max = L, method = "ols")
    x[n, 1] <- est$x.intercept
    x[n, 2:(L + 1)] <- est$ar

#    if (method == "ols") {
#      est<-ar(x[(1+(n-1)*window_size):min(n*window_size,N)],aic=FALSE,order.max = L,method="ols")
#      x_transformed[n,1]<-est$x.intercept
#      x_transformed[n,2:(L+1)]<-est$ar
#    }
#    if (method == "mle") {
#      est<-ar(x[(1+(n-1)*window_size):min(n*window_size,N)],aic=FALSE,order.max = L,method="mle")
#      x_transformed[n,1]<-est$x.mean
#      x_transformed[n,2:(L+1)]<-est$ar
#    }
#    if ((method == "yule-walker") || (method == "yw")) {
#      est<-ar(x[(1+(n-1)*window_size):min(n*window_size,N)],aic=FALSE,order.max = L,method="yule-walker")
#      x_transformed[n,1]<-est$x.mean
#      x_transformed[n,2:(L+1)]<-est$ar
#    }
  }
  return(x)
}
