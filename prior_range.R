#test
a=matrix(rnorm(40,mean=-1,sd=1),nrow=20,ncol=2)
b=matrix(rnorm(120,mean=0,sd=1),nrow=60,ncol=2)
c=matrix(rnorm(40,mean=1,sd=1),nrow=20,ncol=2)
x=rbind(a,b,c)
l1=c(15,25)
l2=c(75,100)
prior_range_x=list(l1,l2)
PriorRange(x,prior_range_x=list(l1,l2),num_init = 1)
#test
#K changepoints 
PriorRangeOrderKmeans=function(x,prior_range_x,num_init=sqrt(dim(x)[1])) {
  if (class(x) != "matrix") {
    stop("Dataset must be matrix form!")
  }
  N=dim(x)[1] # number of observations
  D=dim(x)[2] # dimension of each observation
  #There are K segments 
  K=length(prior_range_x)+1
  # Change points number error handling
  #test
  #print(N)
  #print(K)
  #test
  if (N < K) {
    stop("Change point number too large or Input dimension error!")
  }
  if (N == K) {
    k_changepints=list(wgss=0,changepoints=as.numeric(seq(K)))
    return(k_changepints)
  }
  # randomize initial change points several times to avoid local optima
  best_wgss=Inf
  for (j in 1:num_init) {
    #test
    #cat("j=",j,"\n")
    #test

    # Special case: only have one change point
    if (K==1) {
      changePoints = N
      num_each=N
      wgss=var(x)*(N-1)
      if (N==1) {
        wgss=matrix(0,nrow=1,ncol=D)
      }
    } else {
      # store the within segment sum of squared distances to the segment mean (wgss)
      # in each dimension in each segment
      num_each=matrix(0,nrow=K,ncol=1)
      wgss_each=matrix(0,nrow=K,ncol=D)
      mean_each=matrix(0,nrow=K,ncol=D)

      # initialize change points
      changePoints=numeric(K)
      changePoints[K]=N
      for (i in 1:(K-1)) {
        changePoints[i]=floor(prior_range_x[[i]][1]+(prior_range_x[[i]][2]-prior_range_x[[i]][1])*runif(1))
      }

      # initialize for each segment the number of points, within segment sum of squares, and mean
      num_each[1] = changePoints[1]
      wgss_each[1,] = apply(matrix(x[1:changePoints[1],],ncol=D),2,var) * (num_each[1]-1)
      if (num_each[1]==1) {
        wgss_each[1,]=matrix(0,nrow=1,ncol=D)
      }
      mean_each[1,] = apply(matrix(x[1:changePoints[1],],ncol=D),2,mean)

      for (i in 2:K) {
        num_each[i] = changePoints[i] - changePoints[i-1]
        wgss_each[i,] = apply(matrix(x[(changePoints[i-1]+1):changePoints[i],],ncol=D),2,var) * (num_each[i]-1)
        # special case: one segment only contains one number, avoid get NA for variance
        if (num_each[i]==1) {
          wgss_each[i,]=matrix(0,nrow=1,ncol=D)
        }
        mean_each[i,] = apply(matrix(x[(changePoints[i-1]+1):changePoints[i],],ncol=D),2,mean)

      }
      # scan the middle K-1 change points
      # suppose that we are at the crossing of segments i and i+1
      for (i in 1:(K-1)) {
        #test
        #i=2
        #cat("i=",i,"\n")
        #cat("num_each=",num_each,"\n")
        #cat("num_each[i]=",num_each[i],"\n")
        #cat("changePoints=",changePoints,"\n")
        #test
        
        # consider if we should move the last part of segment i
        best_gain_sum=-Inf
        if ((changePoints[i]-prior_range_x[[i]][1])>0) {

          # scan all possible part that can be transformed form segment i to i+1
          for (ell in 1:(changePoints[i]-prior_range_x[[i]][1])) {
            #test
            #cat("ell=",ell,"\n")
            #test

            mean_candidatePart=apply(matrix(x[(changePoints[i]-ell+1):changePoints[i],],ncol=D),2,mean)
            # the descrease in wgss of segment i
            decrease=ell*num_each[i]/(num_each[i]-ell)*(mean_each[i,]-mean_candidatePart)^2
            # the increase in wgss of segment i+1
            increase=ell*num_each[i+1]/(num_each[i+1]+ell)*(mean_each[i+1,]-mean_candidatePart)^2
            # store the best candidate part than can be transformed from segment i to i+1
            if ( sum(decrease) - sum(increase) > best_gain_sum) {
              #test
              #cat("decrease=",sum(decrease),"\n")
              #cat("increase=",sum(increase),"\n")
              #test
              best_gain = decrease - increase
              best_gain_sum = sum(best_gain)
              best_ell = ell
              best_candidatePart = matrix(x[(changePoints[i]-ell+1):changePoints[i],],ncol=D)
              #test
              #cat("best_candidatePart=",best_candidatePart,"\n")
              #cat("dim_candidate=",dim(best_candidatePart),"\n")
              #cat("best_ell =",best_ell,"\n")
              #test
              best_decrease = decrease
              best_increase = increase
            
            }
          }
        }
        # If total wgss decrease, move the best candidate part from segment i to i+1,
        # and get new segment i and i+1, also update the corresponding mean, wgss adn change point.
        # If not, consider if we should move the first part of segment i+1
        if  (best_gain_sum > 0) {
          #test
          #cat("left to right","\n")
          #test
          best_mean_candidatePart = apply(best_candidatePart,2,mean)
          mean_each[i,] = (num_each[i]*mean_each[i,]-best_ell*best_mean_candidatePart)/(num_each[i]-best_ell)
          mean_each[i+1,] = (num_each[i+1]*mean_each[i+1,]+best_ell*best_mean_candidatePart)/(num_each[i+1]+best_ell)
          changePoints[i] = changePoints[i] - best_ell
          num_each[i] = num_each[i] - best_ell
          num_each[i+1] = num_each[i+1] + best_ell
          wgss_part = apply((best_candidatePart - matrix(best_mean_candidatePart,nrow=best_ell,ncol=D,byrow=TRUE))^2,2,sum)
          wgss_each[i,] = wgss_each[i,] - best_decrease - wgss_part
          wgss_each[i+1,] = wgss_each[i+1,] + best_increase + wgss_part
          #test
          #cat("changePoints[i]",changePoints[i],"\n")
          #cat("num_each[i] =",num_each[i],"\n")
          #cat("num_each[i+1] =",num_each[i+1],"\n")
          #test
        } else {
          # consider if we should move the first part of segment i+1
          best_gain_sum=-Inf
          #test
          #cat("num_each[i+1]=",num_each[i+1],"\n")
          #test
          if ((prior_range_x[[i]][2]-changePoints[i])>0) {
            for (ell in 1:(prior_range_x[[i]][2]-changePoints[i])) {
              #test
              #ell=24
              #cat("ell=",ell,"\n")
              #test
              if (ell<num_each[i+1]) {
                mean_candidatePart=apply(matrix(x[(changePoints[i]+1):(changePoints[i]+ell),],ncol=D),2,mean)
  
                decrease=ell*num_each[i+1]/(num_each[i+1]-ell)*(mean_each[i+1,]-mean_candidatePart)^2
                increase=ell*num_each[i]/(num_each[i]+ell)*(mean_each[i,]-mean_candidatePart)^2
                if ( sum(decrease) - sum(increase) > best_gain_sum) {
                  #test
                  #cat("decrease=",sum(decrease),"\n")
                  #cat("increase=",sum(increase),"\n")
                  #test
                  best_gain = decrease - increase
                  best_gain_sum = sum(best_gain)
                  best_ell = ell
                  best_candidatePart = matrix(x[(changePoints[i]+1):(changePoints[i]+ell),],ncol=D)
                  #test
                  #cat("best_candidatePart=",best_candidatePart,"\n")
                  #cat("dim_candidate=",dim(best_candidatePart),"\n")
                  #cat("best_ell =",best_ell,"\n")
                  #test
                  best_decrease = decrease
                  best_increase = increase
                }
              }
            }
          }
          if  (best_gain_sum > 0) {
            #test
            #cat("right to left","\n")
            #test
            best_mean_candidatePart = apply(best_candidatePart,2,mean)
            mean_each[i+1,] = (num_each[i+1]*mean_each[i+1,]-best_ell*best_mean_candidatePart)/(num_each[i+1]-best_ell)
            mean_each[i,] = (num_each[i]*mean_each[i,]+best_ell*best_mean_candidatePart)/(num_each[i]+best_ell)
            changePoints[i] = changePoints[i] + best_ell
            num_each[i] = num_each[i] + best_ell
            num_each[i+1] = num_each[i+1] - best_ell
            wgss_part = apply((best_candidatePart - matrix(best_mean_candidatePart,nrow=best_ell,ncol=D,byrow=TRUE))^2,2,sum)
            wgss_each[i,] = wgss_each[i,] + best_decrease + wgss_part
            wgss_each[i+1,] = wgss_each[i+1,] - best_increase - wgss_part
            #test
            #cat("changePoints[i]",changePoints[i],"\n")
            #cat("num_each[i] =",num_each[i],"\n")
            #cat("num_each[i+1] =",num_each[i+1],"\n")
            #test
          }
        }
      }
      # get the wgss of all segments with original dimension
      wgss=colSums(wgss_each)
    }
    #get the total wgss of all dimensions
    wgss_sum=sum(wgss)
    # store the smallest total wgss among several initializations.
    if (best_wgss>wgss_sum) {
      best_wgss=wgss_sum
      best_changepoints=changePoints
    }
    #print(wgss_sum)
    #print(changePoints)
  }

  k_changepints=list(num_changepoints=K,changepoints=best_changepoints)
  return(k_changepints)
}

