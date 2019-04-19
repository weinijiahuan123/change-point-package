#test
set.seed(5701)
N = 1000
N1 = floor(0.1*N)
N2 = floor(0.3*N)
L = 2
a1 = c(0.8, -0.3); c1 = 0
a2 = c(-0.5, 0.1); c2 = 0
a3 = c(0.5, -0.5); c3 = 0
x = rep(0,N)
x[1:L] = rnorm(L)

for (n in (L+1):N){
    if (n <= N1) {
        x[n] = x[(n-1):(n-L)] %*% a1 + c1 + rnorm(1)
    }
    else if (n <= (N1+N2)) {
        x[n] = x[(n-1):(n-L)] %*% a2 + c2 + rnorm(1)
    }
    else {
        x[n] = x[(n-1):(n-L)] %*% a3 + c3 + rnorm(1)
    }
}
window_size=10
method="ols"
GetMle(x,window_size=10,L=2,method="ols")
GetMle=function(x,window_size,L,method="ols") {
    # Transform N dependent data into approximated independent data (N/window_size) x (L+1).
    # Computes the estimated coefficients of each window of original data.
    #
    # Args:
    #   x: The original data to find change points.
    #   window_size: The number of observations each window contains.
    #   L: Lag order of the dataset. L>=1
    #   method: The method used to estimate coefficients.
    #           Must be one of c("yule-walker", "own-ols", "ols", "mle", "yw").
    #           Default is ols (ordinal least squares).
    #
    # Returns:
    #   x_transformed: The transformed data, which are the estimated coefficients of original data.
    N=length(x)
    n_window = ceiling(N/window_size)
    x_transformed=matrix(0,nrow=n_window,ncol=L+1)
    for (n in 1:n_window) {
        #test
        #test
        #get estimated coefficients including constant
        if (method == "ols") {
            est=ar(x[(1+(n-1)*window_size):min(n*window_size,N)],aic=FALSE,order.max = L,method="ols")
            x_transformed[n,1]=est$x.intercept
            x_transformed[n,2:(L+1)]=est$ar
        }
        if (method == "mle") {
            est=ar(x[(1+(n-1)*window_size):min(n*window_size,N)],aic=FALSE,order.max = L,method="mle")
            x_transformed[n,1]=est$x.mean
            x_transformed[n,2:(L+1)]=est$ar
        }
        if (method == "yule-walker") {
            est=ar(x[(1+(n-1)*window_size):min(n*window_size,N)],aic=FALSE,order.max = L,method="yule-walker")
            x_transformed[n,1]=est$x.mean
            x_transformed[n,2:(L+1)]=est$ar
        }
        if (method == "own-ols") {

            est=EstimateAr(x,1+(n-1)*window_size,min(n*window_size,N),L)
            #transform original data to transformed data which is the estimated coefficients
            if (n==1) {
                x_transformed=t(est$C)
            } else {
                x_transformed=rbind(x_transformed,t(est$C))
            }
        }
    }
    return(x_transformed)
}




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

test.yw=ar(x[11:20],aic=FALSE,order.max = 2, method="yw")
test=ar(x[11:20],aic=FALSE,order.max = 2, method="mle")
test.ols=ar(c(x[11:20]),aic=FALSE,order.max = 2, method="ols")
x[19]*test.yw$ar[1]+x[18]*test.yw$ar[2]+test.yw$x.mean
