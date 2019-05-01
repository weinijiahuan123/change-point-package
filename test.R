#set.seed(5701)
# simulation 1
#test for OrderKmeans
a=matrix(rnorm(40,mean=-1,sd=1),nrow=20,ncol=2)
b=matrix(rnorm(120,mean=0,sd=1),nrow=60,ncol=2)
c=matrix(rnorm(40,mean=1,sd=1),nrow=20,ncol=2)
x=rbind(a,b,c)
OrderKmeans(x,K=2)
#test for ChangePoints
ChangePoints(x,point_max=5,seg_min = 2)

# simulation 2
N = 1000
N1 = floor(0.1*N)
N2 = floor(0.3*N)
a1 = c(0.8, -0.3); c1 = 0
a2 = c(-0.5, 0.1); c2 = 0
a3 = c(0.5, -0.5); c3 = 0
y = rep(0,N)
L=2
y[1:L] = rnorm(L)
for (n in (L+1):N){
  if (n <= N1) {
    y[n] = y[(n-1):(n-L)] %*% a1 + c1 + rnorm(1)
  } else if (n <= (N1+N2)) {
    y[n] = y[(n-1):(n-L)] %*% a2 + c2 + rnorm(1)
  }
  else {
    y[n] = y[(n-1):(n-L)] %*% a3 + c3 + rnorm(1)
  }
}

#test for multi window
MultiWindow(y,window_list=c(100,50,20,10,5),point_max=5,seg_min=2,tolerance=1)
# If set get_mle= GetMleAr,pay attention to the lag L, sometimes set L=2 will get error:
#  In ar.ols(x, aic = aic, order.max = order.max, na.action = na.action,  :
#              model order:  2 singularities in the computation of the projection matrix results are only valid up to model order 1
MultiWindow(y,window_list=c(100,50,20,10,5),point_max=5,seg_min=1,tolerance=1,get_mle = GetMleAr)
MultiWindow(y,window_list=c(100,50,20,10,5),prior_range=list(c(30,200),c(220,400)))
# test for GetMle
GetMle(y,window_size=100)
# test for GetMleAr
GetMleAr(y,window_size=100)

# simulation 3
N = 1000
N1 = floor(0.1*N)
N2 = floor(0.3*N)
a1 = c(0.8, -0.6); c1 = 0.1
a2 = c(-0.5, 0.5); c2 = 0.2
a3 = c(0.5, -0.7); c3 = 0.3
y = rep(0,N)
L=2
y[1:L] = rnorm(L)
for (n in (L+1):N){
  if (n <= N1) {
    y[n] = y[(n-1):(n-L)] %*% a1 + c1 + rnorm(1)
  } else if (n <= (N1+N2)) {
    y[n] = y[(n-1):(n-L)] %*% a2 + c2 + rnorm(1)
  }
  else {
    y[n] = y[(n-1):(n-L)] %*% a3 + c3 + rnorm(1)
  }
}

#test for multi window
MultiWindow(y,window_list=c(100,50,20,10,5),point_max=5,seg_min=2,tolerance=1)
# If set get_mle= GetMleAr,pay attention to the lag L, sometimes set L=2 will get error:
#  In ar.ols(x, aic = aic, order.max = order.max, na.action = na.action,  :
#              model order:  2 singularities in the computation of the projection matrix results are only valid up to model order 1
MultiWindow(y,window_list=c(100,50,20,10,5),point_max=5,seg_min=1,tolerance=1,get_mle = GetMleAr)

# Simulation 4
N <- 1000
X <- gen_EFdata(N)
res <- MultiWindow(X,window_list=c(100,80,50,30, 20),get_mle=GetHle,point_max=2,seg_min=1,tolerance=1)
n_cp <- res$n_peak_range
loc_cp <- rep(0, n_cp)
for (j in 1:n_cp){
  loc_cp[j] <- mean(res$peak_range[[j]])
}
loc_cp

#install.packages('MHadaptive')
library('MHadaptive')

gen_EFdata <- function(N){
  q <- 4

  Ns <- c(floor(0.2*N), floor(0.3*N), N-floor(0.2*N)-floor(0.3*N))
  thetas <- c(0.2,0.6,1)
  X <- c()
  for (k in 1:3){
    loglik <- function( x ){
      res <- - thetas[k] * abs(x)^q
    }
    mh <- Metro_Hastings(loglik, rgamma(n = 1, shape = 3, rate = 1), prop_sigma = NULL,
                         par_names = NULL, iterations = 11000, burn_in = 1000,
                         adapt_par = c(100, 20, 0.5, 0.75), quiet = TRUE)
    mh0 <-  mcmc_thin(mh, thin = 20)
    #acf(mh0$trace)
    #plotMH(mh0, correlogram = TRUE)
    X <- c(X, mh0$trace[1:Ns[k]])
  }
  #plot(X)
  X
}


GetHle=function(x,window_size) {
  N=length(x)
  n_window = ceiling(N/window_size)
  x_transformed=rep(0,n_window)
  for (n in 1:n_window) {
    #get estimated coefficients
    xx <- x[(1+(n-1)*window_size):min(n*window_size,N)]
    x_transformed[n] <- q * (q-1) * sum(abs(xx)^(q-2)) / sum(q^2 * abs(xx)^(2*q-2))
  }
  return(x_transformed)
}



