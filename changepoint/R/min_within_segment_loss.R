#' Detect Number and Location of Change Points of Independent Data
#'
#' Detect the number and locations of change points based on minimizing within
#' segment quadratic loss and applying penalized model selection approach
#' with restriction of largest candidate number of change points.
#'
#' The K change points form K segments (1 2 ... change_point(1)) ...
#' (change_point(K-1)+1 ... change_point(K)=N).
#'
#' @references
#' J. Ding, Y. Xiang, L. Shen, and V. Tarokh, \emph{Multiple Change Point Analysis:
#' Fast Implementation and Strong Consistency}. IEEE Transactions on Signal
#' Processing, vol. 65, no. 17, pp. 4495-4510, 2017.
#' @param x The data to find change points.
#' @param point_max The largest candidate number of change points.
#' @param penalty Penalty term. Default is log(N).
#' @param seg_min Minimal segment size, must be positive integer.
#' @param num_init The number of repetition times, in order to avoid local minimal.
#'
#' @return A list with following elements:
#'         M_hat: optimal number of change points.
#'         changepoints_hat: location of change points.lag contains the number and location of change points.
#' @export
#' @examples
#' a=matrix(rnorm(40,mean=-1,sd=1),nrow=20,ncol=2)
#' b=matrix(rnorm(120,mean=0,sd=1),nrow=60,ncol=2)
#' c=matrix(rnorm(40,mean=1,sd=1),nrow=20,ncol=2)
#' x=rbind(a,b,c)
#' ChangePoints(x,point_max=5)
#' ChangePoints(x,point_max=5,penalty="hq")

ChangePoints<-function(x,point_max=4,penalty="bic",seg_min=2,num_init="sqrt"){
  N<-dim(x)[1]
  D<-dim(x)[2]
  #sigma2<-sum(apply(x,2,var))
  #wgss_list=list()
  #wgss_list_penalty=list()
  wgss_llist_sigma=list()
  #wgss_llist_2sigma=list()
  best_wgss_penalty<-Inf
  # Make sure the number of change points is no larger than the number of observations.
  point_max<-min(point_max,N)

  # K = number of segments

  for (K in 1:point_max) {
    #K=2
    res <- OrderKmeans(x,K,num_init=num_init)
    wgss <- res$wgss
    num_each <- res$num_each
    changepoints <- res$changepoints
    #if size of the smallest segment is less than minimal segment size, end the loop
    if (min(num_each)<seg_min) {
      break
    }
    #get the within segment sum of residual plus penalty
    #wgss_list=c(wgss_list,wgss)
    #wgss_penalty_new=wgss+K*penalty
    #wgss_penalty<-wgss/sigma2+K*penalty
    if (penalty == "bic") {
      wgss_penalty <- sum(num_each * log(wgss/num_each)) + D * K * log(N)
    } else if (penalty == "aic") {
      wgss_penalty <- sum(num_each * log(wgss/num_each)) + 2 * D * K
    } else if (penalty == "hq") {
      wgss_penalty <- sum(num_each * log(wgss/num_each)) + D * K * log(log(N))
    }
    #wgss_list_penalty=c(wgss_list_penalty,wgss_penalty_new)
    wgss_llist_sigma=c(wgss_llist_sigma,wgss_penalty)
    #wgss_llist_2sigma=c(wgss_llist_2sigma,wgss/(2*sigma2)+K*penalty)
    #test
    #print(penalty)
    #test

    #if the new wgss plus penalty is smaller than the previous one, store the new value,
    #get the smallest wgss among different changepoint numbers
    if (wgss_penalty<best_wgss_penalty) {
      best_wgss_penalty<-wgss_penalty
      best_changepoints<-changepoints
      best_num<-K
    }
  }
  # test
  #print(wgss_list)
  #print(wgss_list_penalty)
  print(wgss_llist_sigma)
  #print(wgss_llist_2sigma)
  # test
  m_changepoints<-list(num_changepoints=best_num,changepoints=best_changepoints)
  return(m_changepoints)
}

#' Detect Location of Change Points of Independent Data
#'
#' Detect the location of change points based on minimizing within segment quadratic
#' loss with fixed number of change points.
#'
#' The K change points form K segments (1 2 ... change_point(1)) ...
#' (change_point(K-1)+1 ... change_point(K)=N).
#'
#' @references
#' J. Ding, Y. Xiang, L. Shen, and V. Tarokh, \emph{Multiple Change Point Analysis:
#' Fast Implementation and Strong Consistency}. IEEE Transactions on Signal
#' Processing, vol. 65, no. 17, pp. 4495-4510, 2017.
#' @param x The data to find change points with dimension N x D, must be matrix
#' @param K The number of change points.
#' @param num_init The number of repetition times, in order to avoid local minimal.
#'                 Default is squared root of number of observations
#'
#' @return A list contains following elements:
#'         wgss: smallest within segment sum of squared distances to the segment mean.
#'         changepoints: location of optimal change points.
#' @export
#' @examples
#' a=matrix(rnorm(40,mean=-1,sd=1),nrow=20,ncol=2)
#' b=matrix(rnorm(120,mean=0,sd=1),nrow=60,ncol=2)
#' c=matrix(rnorm(40,mean=1,sd=1),nrow=20,ncol=2)
#' x=rbind(a,b,c)
#' OrderKmeans(x,K=3)
#' OrderKmeans(x,K=3,num_init=2*sqrt(dim(x)[1]))
OrderKmeans <- function(x, K=4, num_init="sqrt") {
  if (class(x) != "matrix") {
    stop("Dataset must be matrix form!")
  }
  N<-dim(x)[1] # number of observations
  D<-dim(x)[2] # dimension of each observation
  # There are M segments
  M <- K+1
  # Change points number error handling
  if (N < M) {
    stop("Change point number too large or Input dimension error!")
  } else if (M == 1) {
    stop("No change point!")
  }
  # special case, N == M
  if (N == M) {
    best_wgss_sum <- 0
    best_num_each <- matrix(1,nrow=N,ncol=1)
    best_wgss <- rep(0,N)
    best_changepoints <- as.numeric(seq(N))
  } else {

    # randomize initial change points several times to avoid local optima
    best_wgss_sum<-Inf
    # set the random initialization times
    if (num_init == "sqrt") {
      num_init <- sqrt(dim(x)[1])
    }
    for (j in 1:num_init) {
      #test
      #cat("j=",j)
      #test

      # store the within segment sum of squared distances to the segment mean (wgss)
      # in each dimension in each segment
      num_each<-matrix(0,nrow=M,ncol=1)
      wgss_each<-matrix(0,nrow=M,ncol=D)
      mean_each<-matrix(0,nrow=M,ncol=D)

      # Initialize change points, add N-th observation as the last change point so
      # there are K+1 change points
      changePoints<-floor(1+(N-2)*runif(K))
      changePoints<-c(changePoints,N)
      changePoints<-unique(changePoints)

      # make sure change points are unique
      while (length(changePoints)<M) {
        changePoints<-c(changePoints,floor(1+(N-2)*runif(1)))
        changePoints<-unique(changePoints)
      }
      changePoints <- sort(changePoints)

      # initialize for each segment the number of points, within segment sum of squares, and mean
      num_each[1] <- changePoints[1]
      wgss_each[1,] <- apply(matrix(x[1:changePoints[1],],ncol=D),2,var) * (num_each[1]-1)
      if (num_each[1]==1) {
        wgss_each[1,]<-matrix(0,nrow=1,ncol=D)
      }
      mean_each[1,] <- apply(matrix(x[1:changePoints[1],],ncol=D),2,mean)

      for (i in 2:M) {
        num_each[i] <- changePoints[i] - changePoints[i-1]
        wgss_each[i,] <- apply(matrix(x[(changePoints[i-1]+1):changePoints[i],],ncol=D),2,var) * (num_each[i]-1)
        # special case: one segment only contains one number, avoid get NA for variance
        if (num_each[i]==1) {
          wgss_each[i,]<-matrix(0,nrow=1,ncol=D)
        }
        mean_each[i,] <- apply(matrix(x[(changePoints[i-1]+1):changePoints[i],],ncol=D),2,mean)

      }

      # scan the middle K change points
      # suppose that we are at the crossing of segments i and i+1
      for (i in 1:K) {
        #test
        #i=1
        #cat("i=",i,"\n")
        #cat("num_each=",num_each,"\n")
        #cat("num_each[i]=",num_each[i],"\n")
        #test

        # consider if we should move the last part of segment i
        best_gain_sum <- -Inf
        if (num_each[i]>1) {

          # scan all possible part that can be transformed form segment i to i+1
          for (ell in 1:(num_each[i]-1)) {
            #test
            #ell=12
            #cat("ell=",ell,"\n")
            #test

            mean_candidatePart<-apply(matrix(x[(changePoints[i]-ell+1):changePoints[i],],ncol=D),2,mean)
            # the descrease in wgss of segment i
            decrease<-ell*num_each[i]/(num_each[i]-ell)*(mean_each[i,]-mean_candidatePart)^2
            # the increase in wgss of segment i+1
            increase<-ell*num_each[i+1]/(num_each[i+1]+ell)*(mean_each[i+1,]-mean_candidatePart)^2
            # store the best candidate part than can be transformed from segment i to i+1
            if ( sum(decrease) - sum(increase) > best_gain_sum) {
              #test
              #cat("decrease=",sum(decrease),"\n")
              #cat("increase=",sum(increase),"\n")
              #test
              best_gain <- decrease - increase
              best_gain_sum <- sum(best_gain)
              best_ell <- ell
              best_candidatePart <- matrix(x[(changePoints[i]-ell+1):changePoints[i],],ncol=D)
              #test
              #cat("best_candidatePart=",best_candidatePart,"\n")
              #cat("best_ell =",best_ell,"\n")
              #test
              best_decrease <- decrease
              best_increase <- increase
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
          best_mean_candidatePart <- apply(best_candidatePart,2,mean)
          mean_each[i,] <- (num_each[i]*mean_each[i,]-best_ell*best_mean_candidatePart)/(num_each[i]-best_ell)
          mean_each[i+1,] <- (num_each[i+1]*mean_each[i+1,]+best_ell*best_mean_candidatePart)/(num_each[i+1]+best_ell)
          changePoints[i] <- changePoints[i] - best_ell
          num_each[i] <- num_each[i] - best_ell
          num_each[i+1] <- num_each[i+1] + best_ell
          wgss_part <- apply((best_candidatePart - matrix(best_mean_candidatePart,nrow=best_ell,ncol=D,byrow=TRUE))^2,2,sum)
          wgss_each[i,] <- wgss_each[i,] - best_decrease - wgss_part
          wgss_each[i+1,] <- wgss_each[i+1,] + best_increase + wgss_part
          #if (i == 1) {
          #  wgss_each[i,] <- apply(matrix(x[1:changePoints[i],],ncol=D),2,var) * (num_each[i]-1)
          #} else {
          #  wgss_each[i,] <- apply(matrix(x[(changePoints[i-1]+1):changePoints[i],],ncol=D),2,var) * (num_each[i]-1)
          #}
          #wgss_each[i+1,] <- apply(matrix(x[(changePoints[i]+1):changePoints[i+1],],ncol=D),2,var) * (num_each[i+1]-1)
          #test
          #cat("changePoints[i]",changePoints[i],"\n")
          #cat("num_each[i] =",num_each[i],"\n")
          #cat("num_each[i+1] =",num_each[i+1],"\n")
          #test
        } else {
          # consider if we should move the first part of segment i+1
          best_gain_sum<- -Inf
          #test
          #cat("num_each[i+1]=",num_each[i+1],"\n")
          #test
          if (num_each[i+1]>1) {
            for (ell in 1:(num_each[i+1]-1)) {
              #test
              #cat("ell=",ell,"\n")
              #test
              mean_candidatePart<-apply(matrix(x[(changePoints[i]+1):(changePoints[i]+ell),],ncol=D),2,mean)

              decrease<-ell*num_each[i+1]/(num_each[i+1]-ell)*(mean_each[i+1,]-mean_candidatePart)^2
              increase<-ell*num_each[i]/(num_each[i]+ell)*(mean_each[i,]-mean_candidatePart)^2
              if ( sum(decrease) - sum(increase) > best_gain_sum) {
                #test
                #cat("decrease=",sum(decrease),"\n")
                #cat("increase=",sum(increase),"\n")
                #test
                best_gain <- decrease - increase
                best_gain_sum <- sum(best_gain)
                best_ell <- ell
                best_candidatePart <- matrix(x[(changePoints[i]+1):(changePoints[i]+ell),],ncol=D)
                #test
                #cat("best_candidatePart=",best_candidatePart,"\n")
                #cat("best_ell =",best_ell,"\n")
                #test
                best_decrease <- decrease
                best_increase <- increase
              }
            }
          }
          if  (best_gain_sum > 0) {
            #test
            #cat("right to left","\n")
            #test
            best_mean_candidatePart = apply(best_candidatePart,2,mean)
            mean_each[i+1,] <- (num_each[i+1]*mean_each[i+1,]-best_ell*best_mean_candidatePart)/(num_each[i+1]-best_ell)
            mean_each[i,] <- (num_each[i]*mean_each[i,]+best_ell*best_mean_candidatePart)/(num_each[i]+best_ell)
            changePoints[i] <- changePoints[i] + best_ell
            num_each[i] <- num_each[i] + best_ell
            num_each[i+1] <- num_each[i+1] - best_ell
            wgss_part <- apply((best_candidatePart - matrix(best_mean_candidatePart,nrow=best_ell,ncol=D,byrow=TRUE))^2,2,sum)
            wgss_each[i,] <- wgss_each[i,] + best_decrease + wgss_part
            wgss_each[i+1,] <- wgss_each[i+1,] - best_increase - wgss_part
            #test
            #cat("changePoints[i]",changePoints[i],"\n")
            #cat("num_each[i] =",num_each[i],"\n")
            #cat("num_each[i+1] =",num_each[i+1],"\n")
            #test
          }
        }
        # get the total wgss of each segment
        wgss<-rowSums(wgss_each)
      }
      #get the total wgss of all segments
      wgss_sum<-sum(wgss)
      # store the smallest total wgss among several initializations.
      if (best_wgss_sum>wgss_sum) {
        best_wgss_sum<-wgss_sum
        best_wgss <- wgss
        #test
        #best_wgss_test<-wgss_each
        #test
        best_num_each <- num_each
        best_changepoints<-changePoints
      }
      #print(wgss_sum)
      #print(changePoints)
    }
  }
  # Delete the last observation from change points
  best_changepoints=best_changepoints[-M]

  k_changepints<-list(wgss_sum=best_wgss_sum,num_each=best_num_each,wgss=best_wgss,changepoints=best_changepoints)
  return(k_changepints)
}

