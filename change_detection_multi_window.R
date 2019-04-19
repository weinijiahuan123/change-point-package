MultiWindow=function(x,window_list=c(100,50,20,10,5),point_max=5,L=2,penalty=log(N),seg_min=1,num_init=sqrt(N),tolerance=1, method="ols") {
    # Use a sequence of window sizes to capture ranges of change points.
    # Transform N dependent data into approximated independent data (N/window_size) x (L+1).
    # Computes the estimated coefficients of each window of original data.
    #
    # Args:
    #   x: The original data to find change points. Must be one dimensional data
    #   window_list: The list of window sizes, must be in form c(100,50,20,10,5),
    #                in descending order and each window_size > 2L.
    #   point_max: The largest candidate number of change points.
    #   L: Lag order of the dataset. L>=1
    #   penalty: Penalty term
    #   seg_min: Minimal segment size, must be positive integer.
    #   num_init: The number of repetition times, in order to avoid local minimal.
    #   tolerance: The tolerance level , the selected narrow ranges are with score at least S-tolerance
    #   method: Character string giving the method used to estimate coefficients. Default is ols (ordinal least squares).
    #
    # Returns:
    #   n_peak_range: The number of peak ranges
    #   peak_ranges: The location of peak ranges
    len=length(x)
    n_window_type = length(window_list)
    #initialize score matrix
    score=matrix(0,nrow=len,ncol=n_window_type)
    for (r in 1:n_window_type) {
        #test
        #r=5
        #test
        window_size = window_list[r]
        n_window = ceiling(len/window_size)
        # Get transformed approximated independent data
        x_transformed=GetMle(x,window_size=window_size,L=L,method="ols")
        # if the data is a list, transform it into matrix
        if (class(x_transformed) != "matrix") {
            x_transformed=as.matrix(x_transformed)
        }
        # Get the change points of transformed data
        changePoints=ChangePoints(x_transformed,point_max=point_max,penalty=log(log(N)),seg_min=1,num_init=sqrt(N))$changepoints
        # Map the change points of transformed data to original data and get score the change points.
        # don't score the last number of change points, which is the last number of transformed data
        if (length(changePoints) == 1) {
            if (r!=1) {
                score[1:len,r]=score[1:len,r-1]
            }
        } else {
            if (r==1){
                for (k in 1:(length(changePoints)-1)) {
                    score[(1+(changePoints[k]-1)*window_size):min((changePoints[k]+1)*window_size,len),r]=score[(1+(changePoints[k]-1)*window_size):min((changePoints[k]+1)*window_size,len),r]+1
                }
            } else {
                score[1:len,r]=score[1:len,r-1]
                for (k in 1:(length(changePoints)-1)) {
                    score[(1+(changePoints[k]-1)*window_size):min((changePoints[k]+1)*window_size,len),r]=score[(1+(changePoints[k]-1)*window_size):min((changePoints[k]+1)*window_size,len),(r-1)]+1
                }
            }
        }
    }
    peakranges=PeakRange(score=score,tolerance=tolerance,point_max=point_max)
    return(peakranges)
    #return(score)
}
peakrange=MultiWindow(x,window_list=c(100,50,20,10,5),point_max=3,L=2,penalty=log(log(N)),
                     seg_min=1,num_init=sqrt(N),tolerance=1, method="ols")


window_list=c(100,50,20,10,5)
point_max=4
penalty=log(N)
L=2
num_init=sqrt(N)
seg_min=2
tolerance=1
changePoints=c(1,5,10)
set.seed(5701)
a=rnorm(1,mean=-1,sd=1)
b=rnorm(3,mean=0,sd=1)
c=rnorm(1,mean=1,sd=1)
#X=matrix(c(a,b,c))
x=c(a,b,c)
x_transformed=x
    len=length(x)
    n_window_type = length(window_list)
    #initialize score matrix
    score=matrix(0,nrow=len,ncol=n_window_type)
    for (r in 1:n_window_type) {
        #test
        r=5
        #test
        window_size = window_list[r]
        n_window = ceiling(len/window_size)
        # Get transformed approximated independent data
        x_transformed=GetMle(x,window_size=window_size,L=2,method="ols")
        # if the data is a list, transform it into matrix
        if (class(x_transformed) != "matrix") {
            x_transformed=as.matrix(x_transformed)
        }
        # Get the change points of transformed data
        changePoints=ChangePoints(x_transformed,point_max=5,penalty=log(log(N)),seg_min=1,num_init=sqrt(N))$changepoints
        # Map the change points of transformed data to original data and get score the change points.
        # don't score the last number of change points, which is the last number of transformed data
        if (length(changePoints) == 1) {
            if (r!=1) {
                score[1:len,r]=score[1:len,r-1]
            }
        } else {
            if (r==1){
                for (k in 1:(length(changePoints)-1)) {
                    score[(1+(changePoints[k]-1)*window_size):min((changePoints[k]+1)*window_size,len),r]=score[(1+(changePoints[k]-1)*window_size):min((changePoints[k]+1)*window_size,len),r]+1
                }
            } else {
                score[1:len,r]=score[1:len,r-1]
                for (k in 1:(length(changePoints)-1)) {
                    score[(1+(changePoints[k]-1)*window_size):min((changePoints[k]+1)*window_size,len),r]=score[(1+(changePoints[k]-1)*window_size):min((changePoints[k]+1)*window_size,len),(r-1)]+1
                }
            }
        }
    }
    #peakranges=PeakRange(score,tolerance=1,M_max=4)
    #return(peakranges)
    return(score)
