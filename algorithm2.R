set.seed(5701)
a=matrix(rnorm(40,mean=-1,sd=1),nrow=20,ncol=2)
b=matrix(rnorm(120,mean=0,sd=1),nrow=60,ncol=2)
c=matrix(rnorm(40,mean=1,sd=1),nrow=20,ncol=2)
#X=matrix(c(a,b,c))
X=rbind(a,b,c)
N=dim(X)[1];D=dim(X)[2]
M_max=4
num_min=2
wgss_list=list()
wgss_penalty_old=Inf
for (K in 1:M_max) {
if (N <= K) {
    stop("Input dimension error!")
}
#test
#penalty=log(N)
penalty=0
#test
if (K==1) {
    changePoints = N
    num_each=N
    wgss=var(X)*(N-1)
    if (N==1) {
        wgss=matrix(0,nrow=1,ncol=D)
    }
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
    wgss_each[1,] = apply(matrix(X[1:changePoints[1],],ncol=D),2,var) * (num_each[1]-1)
    if (num_each[1]==1) {
        wgss_each[1,]=matrix(0,nrow=1,ncol=D)
    }
    mean_each[1,] = apply(matrix(X[1:changePoints[1],],ncol=D),2,mean)
    
    for (i in 2:K) {
        num_each[i] = changePoints[i] - changePoints[i-1]
        # segment 1 only has one number with var NA
        wgss_each[i,] = apply(matrix(X[(changePoints[i-1]+1):changePoints[i],],ncol=D),2,var) * (num_each[i]-1)
        if (num_each[i]==1) {
            wgss_each[i,]=matrix(0,nrow=1,ncol=D)
        }
        mean_each[i,] = apply(matrix(X[(changePoints[i-1]+1):changePoints[i],],ncol=D),2,mean) 
        
    }
    
    iter=0; move=1; maxIter = N
    while ((move==1) & (iter < maxIter)) {
        move=0
        iter=iter+1
        # test
        #i=1;
        cat("iter=",iter,"\n")
        #test
        for (i in 1:(K-1)) {
            #test 
            cat("i=",i,"\n")
            cat("num_each=",num_each,"\n")
            cat("num_each[i]=",num_each[i],"\n")
            #test
            best_gain_sum=-Inf
            if (num_each[i]>1) {
                
                for (ell in 1:(num_each[i]-1)) { 
                    #test
                    cat("ell=",ell,"\n")
                    #test
                    
                    mean_candidatePart=apply(matrix(X[(changePoints[i]-ell+1):changePoints[i],],ncol=D),2,mean)
                    
                    decrease=ell*num_each[i]/(num_each[i]-ell)*(mean_each[i,]-mean_candidatePart)^2 
                    increase=ell*num_each[i+1]/(num_each[i+1]+ell)*(mean_each[i+1,]-mean_candidatePart)^2
             
                    if ( sum(decrease) - sum(increase) > best_gain_sum) { 
                        #test
                        cat("decrease=",sum(decrease),"\n")
                        cat("increase=",sum(increase),"\n")
                        #test
                        best_gain = decrease - increase
                        best_gain_sum = sum(best_gain)
                        best_ell = ell
                        best_candidatePart = matrix(X[(changePoints[i]-ell+1):changePoints[i],],ncol=D) 
                        #test
                        cat("best_candidatePart=",best_candidatePart,"\n")
                        cat("best_ell =",best_ell,"\n")
                        #test
                        best_decrease = decrease
                        best_increase = increase
                    }
                }
            }
            if  (best_gain_sum > 0) {
                #test
                cat("left to right","\n")
                #test
                move = 1
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
                cat("changePoints[i]",changePoints[i],"\n")
                cat("num_each[i] =",num_each[i],"\n")
                cat("num_each[i+1] =",num_each[i+1],"\n")
                #test
            } else {
                best_gain_sum=-Inf
                #test 
                cat("num_each[i+1]=",num_each[i+1],"\n")
                #test
                if (num_each[i+1]>1) {
                    for (ell in 1:(num_each[i+1]-1)) { 
                        #test
                        cat("ell=",ell,"\n")
                        #test
                        mean_candidatePart=apply(matrix(X[(changePoints[i]+1):(changePoints[i]+ell),],ncol=D),2,mean)
                        
                        decrease=ell*num_each[i+1]/(num_each[i+1]-ell)*(mean_each[i+1,]-mean_candidatePart)^2 
                        increase=ell*num_each[i]/(num_each[i]+ell)*(mean_each[i,]-mean_candidatePart)^2 
                        if ( sum(decrease) - sum(increase) > best_gain_sum) { 
                            #test
                            cat("decrease=",sum(decrease),"\n")
                            cat("increase=",sum(increase),"\n")
                            #test
                            best_gain = decrease - increase
                            best_gain_sum = sum(best_gain)
                            best_ell = ell
                            best_candidatePart = matrix(X[(changePoints[i]+1):(changePoints[i]+ell),],ncol=D) 
                            #test
                            cat("best_candidatePart=",best_candidatePart,"\n")
                            cat("best_ell =",best_ell,"\n")
                            #test
                            best_decrease = decrease
                            best_increase = increase
                        }
                    }
                }
                if  (best_gain_sum > 0) {
                    #test
                    cat("right to left","\n")
                    #test
                    move = 1
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
                    cat("changePoints[i]",changePoints[i],"\n")
                    cat("num_each[i] =",num_each[i],"\n")
                    cat("num_each[i+1] =",num_each[i+1],"\n")
                    #test
                }
            }
        }
        #test
        cat("move=",move,"\n")
        #test
    }
    wgss=colSums(wgss_each)
}
if (min(num_each)<num_min) {
    break
}
wgss_sum_new=sum(wgss)
wgss_list=c(wgss_list,wgss_sum_new)
wgss_penalty_new=wgss_sum_new+K*penalty
if (wgss_penalty_new<wgss_penalty_old) {
    wgss_penalty_old=wgss_penalty_new
    changePoints_optimal=changePoints
}
}

traceback()
