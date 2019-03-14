#test
set.seed(5701)
N = 2000
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
windowSizes = c(300, 250, 200, 150, 100, 50)
MultiWindow(x,windowSizes=c(300, 250, 200, 150, 100, 50),M_max=4,penalty=log(N),num_min=1,n_init=3,tolerance=1)
#test

MultiWindow=function(x,windowSizes=c(300, 250, 200, 150, 100, 50),M_max=4,penalty=log(N),num_min=1,n_init=3,tolerance=1) {
    N=length(x)
    nWindowType = length(windowSizes)
    score=matrix(0,nrow=N,ncol=nWindowType)
    for (r in 1:nWindowType) {
        #test
        #r=4
        #test
        windowSize = windowSizes[r]
        nWindow = ceiling(N/windowSize)
        for (n in 1:nWindow) {
            #test
            #n=1
            #test
            est=estimate_ar(x,1+(n-1)*windowSize,min(n*windowSize,N),L)
            if (n==1) {
                x_transformed=t(est$C)
            } else {
            x_transformed=rbind(x_transformed,t(est$C))
            }
        }
        changePoints=OrderKmeans(x_transformed,M_max=4,penalty=log(N),num_min=1,n_init=3)$changepoints_hat
        if (r==1){
            for (k in 1:(length(changePoints))) {
                score[(1+(changePoints[k]-1)*windowSize):min((changePoints[k]+1)*windowSize,N),r]=score[(1+(changePoints[k]-1)*windowSize):min((changePoints[k]+1)*windowSize,N),r]+1
            }
        } else {
            for (k in 1:(length(changePoints))) {
                score[1:N,r]=score[1:N,r-1]
                score[(1+(changePoints[k]-1)*windowSize):min((changePoints[k]+1)*windowSize,N),r]=score[(1+(changePoints[k]-1)*windowSize):min((changePoints[k]+1)*windowSize,N),(r-1)]+1
            }
        }
    }
    peakranges=PeakRange(score,tolerance=1,M_max=4)
    return(peakranges)
}

#test
y=x
T1=1+(k-1)*windowSize
T2=min(k*windowSize,N)
L=2
#test
estimate_ar=function(y,T1,T2,L){
    if (T1>(T2-L)) {
        warning("Error in estimate_ar")
    }
    if (T1<=L) {
        T1=L+1
    }
    Y=matrix(0,nrow=L+1,ncol=T2-T1+1)
    Y[1,]=1
    for (k in 1:L) {
        Y[k+1,]=y[(T1-k):(T2-k)]
    }
    A=Y%*%t(Y)
    B=Y%*%y[T1:T2]
    C=solve(A)%*%B
    e=y[T1:T2]-t(Y)%*%C
    sigma2=sum(e^2)/(T2-T1+1)
    est_coef=list(C=C,sigma2=sigma2)
    return(est_coef)
}
#test
k=1
#test
est=estimate_ar(x,1+(k-1)*windowSize,min(k*windowSize,N),L)

