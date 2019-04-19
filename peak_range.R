PeakRange=function(score,tolerance=1,point_max=4) {
    # Select the narrow peak ranges. For each column(window type), find the union of all the peak
    # ranges whose associated scores are no less than S - tolerance, where S is highest score, then
    #
    # Args:
    #   score: The score data to peak ranges.
    #   tolerance: The tolerance level , the selected narrow ranges are with score at least S-tolerance
    #   point_max: The largest candidate number of change points.
    #
    # Returns:
    #   n_peak_range: The number of peak ranges
    #   peak_ranges: The location of peak ranges
    N=dim(score)[1]
    R=dim(score)[2]
    # compute the hightest score
    S=max(score)
    for (r in R:1) {
        J=list()
        # store number of unions
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
                # update the number of unions when meet new union
                if (score[i,r]!=score[i-1,r]) {
                    num=num+1
                    J[[num]]=i
                } else {
                    J[[num]]=c(J[[num]],i)
                }
            }

        }
        # when the number of unions meet the requirement of largest number of change points, end the loop
        if (num<=point_max) {
            break
        }
    }
    optimal=list(n_peak_range=num,peak_ranges=J)
    return(optimal)
}
#test
PeakRange(score,tolerance=1,point_max=5)
