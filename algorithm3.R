score=matrix(c(0,0,1,2,0,1,2,2,1,2,3,4,1,2,3,3,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,2,0,1,2,3,0,1,2,3,0,0,0,0,0,0,0,0),ncol=4,byrow=TRUE)
tolerance=2

N=dim(score)[1];R=dim(score)[2]
S=max(score[,-1])
M_max=3
PeakRange=function(score,tolerance=1,M_max=4) {
    N=dim(score)[1];R=dim(score)[2]
    S=max(score[,-1])
    for (r in R:1) {
        J=list()
        num=0
        #test
        #r=4
        #test
        if (score[1,r]>=S-tolerance) {
            num=num+1
            J[[num]]=1
        }
        for (i in 2:N) {
            if (score[i,r]>=S-tolerance) {
                if (score[i,r]!=score[i-1,r]) {
                    num=num+1
                    J[[num]]=i
                } else {
                    J[[num]]=c(J[[num]],i)
                }
            }
     
        }
        if (num<=M_max) {
            break
        }
    }
    optimal=list(best_nRange=num,best_Ranges=J,best_nWindowType=r)
    return(optimal)
}
#output
M_optimal=num
cp=J
r_optimal=r
#x_optimal=list()
#for (i in 1:num) {
#    print(score[J[[i]],r])
#}

