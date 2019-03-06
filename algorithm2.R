set.seed(5701)
a=rnorm(20,mean=-1,sd=1)
b=rnorm(60,mean=0,sd=1)
c=rnorm(20,mean=1,sd=1)
X=matrix(c(a,b,c))
N=dim(X)[1];D=dim(X)[2]
K=5
if (N <= K) {
  warning("Input dimension error!")
} 
if (K==1) {
  changePoints = N
  wgss=var(x)*(N-1)
} else {
  num_each=matrix(0,nrow=K,ncol=1)
  wgss_each=matrix(0,nrow=K,ncol=D)
  mean_each=matrix(0,nrow=K,ncol=D)
  changePoints=floor(1+(N-2)*runif(K-1))
  changePoints=c(changePoints,N)
  changePoints=unique(changePoints)
  
  while (length(changePoints)<K) {
    changePoints=c(changePoints,floor(1+(N-2)*runif(1)))
    changePoints=unique(changePoints)
  }
  
  changePoints = sort(changePoints)
  
  num_each[1] = changePoints[1] 
  wgss_each[1,] = apply(matrix(X[1:changePoints[1],]),2,var) * (num_each[1]-1)
  mean_each[1,] = apply(matrix(X[1:changePoints[1],]),2,mean)
  
  for (i in 2:K) {
  num_each[i] = changePoints[i] - changePoints[i-1]
  # segment 4 only has one number with var NA
  wgss_each[i,] = apply(matrix(X[(changePoints[i-1]+1):changePoints[i],]),2,var) * (num_each[i]-1)
  mean_each[i,] = apply(matrix(X[(changePoints[i-1]+1):changePoints[i],]),2,mean) 
  }
  
  iter=0; move=1; maxIter = N
  while ((move==1) & (iter < maxIter)) {
    move=0
    iter=iter+1
    # test
    #i=1
    #test
    for (i in 1:K-1) {
      best_gain_sum=-Inf
      for (ell in 1:num_each[i]-1) { 
        #test
        #ell=2
        #test
        mean_candidatePart=apply(matrix(X[(changePoints[i]-ell+1):changePoints[i],]),2,mean)
      
        decrease=ell*num_each[i]/(num_each[i]-ell)*(mean_each[i,]-mean_candidatePart)^2 
        increase=ell*num_each[i+1]/(num_each[i+1]+ell)*(mean_each[i+1,]-mean_candidatePart)^2 
        if ( sum(decrease) - sum(increase) > best_gain_sum) { 
          best_gain = decrease - increase
          best_gain_sum = sum(best_gain)
          best_ell = ell
          best_candidatePart = matrix(X[(changePoints[i]-ell+1):changePoints[i],]) 
          best_decrease = decrease
          best_increase = increase
        }
      }
      if  (best_gain_sum > 0) {
        move = 1
        best_mean_candidatePart = apply(best_candidatePart,2,mean)
        mean_each[i,] = (num_each[i]*mean_each[i,]-best_ell*best_mean_candidatePart)/(num_each[i]-best_ell) 
        mean_each[i+1,] = (num_each[i+1]*mean_each[i+1,]+best_ell*best_mean_candidatePart)/(num_each[i+1]+best_ell) 
        changePoints[i] = changePoints[i] - best_ell
        num_each[i] = num_each[i] - best_ell
        num_each[i+1] = num_each[i+1] + best_ell
        wgss_part = apply((best_candidatePart - matrix(best_mean_candidatePart,nrow=best_ell,ncol=D,byrow=TRUE))^2,2,sum)
        wgss_each[i,] = wgss_each[i,] - best_decrease - wgss_part
        wgss_each[i+1,] = wgss_each[i+1] + best_increase + wgss_part  
      } else {
        best_gain_sum=-Inf
        for (ell in 1:num_each[i+1]-1) { 
          mean_candidatePart=apply(matrix(X[(changePoints[i]+1):changePoints[i]+ell,]),2,mean)
          
          decrease=ell*num_each[i+1]/(num_each[i+1]-ell)*(mean_each[i+1,]-mean_candidatePart)^2 
          increase=ell*num_each[i]/(num_each[i]+ell)*(mean_each[i,]-mean_candidatePart)^2 
          if ( sum(decrease) - sum(increase) > best_gain_sum) { 
            best_gain = decrease - increase
            best_gain_sum = sum(best_gain)
            best_ell = ell
            best_candidatePart = matrix(X[(changePoints[i]+1):changePoints[i]+ell,]) 
            best_decrease = decrease
            best_increase = increase
          }
        }
        if  (best_gain_sum > 0) {
          move = 1
          best_mean_candidatePart = apply(best_candidatePart,2,mean)
          mean_each[i+1,] = (num_each[i+1]*mean_each[i+1,]-best_ell*best_mean_candidatePart)/(num_each[i+1]-best_ell) 
          mean_each[i,] = (num_each[i]*mean_each[i,]+best_ell*best_mean_candidatePart)/(num_each[i]+best_ell) 
          changePoints[i] = changePoints[i] + best_ell
          num_each[i] = num_each[i] + best_ell
          num_each[i+1] = num_each[i+1] - best_ell
          wgss_part = apply((best_candidatePart - matrix(best_mean_candidatePart,nrow=best_ell,ncol=D,byrow=TRUE))^2,2,sum)
          wgss_each[i,] = wgss_each[i,] + best_decrease + wgss_part
          wgss_each[i+1,] = wgss_each[i+1] - best_increase - wgss_part  
        }
      }
    }
    wgss=colSums(wgss_each)
  }
}
wgss_sum=sum(wgss)

