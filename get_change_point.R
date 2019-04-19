RangeToPoint=function(x,n,location) {
    # Transform the peak ranges of change points to exact change points.
    #
    # Args:
    #   x: The original data to find change points. Must be one dimensional data
    #   n: The number of peak ranges of change points
    #   location: The location of ranges of change points
    #
    # Returns:
    #   change_point: The location of exact change points
    len=length(x)
    change_point=numeric(n+1)
    # add 1 as the first change point for computation convenience
    change_point[1]=1
    for (i in 1:n) {
        best_loss=Inf
        #test
        #i=1
        #test
        #cat("i=",i,"\n")
        # Get the outer range to find the exact change point.
        if (n==1) {
            outer_range=x
        } else {
            if (i == n) {
                outer_range=x[change_point[i]:len]
            } else {
            outer_range=x[change_point[i]:(location[[i+1]][1]-1)]
            }
        }
        # Fix following ranges and find the exact change point of range i with smallest quadratic loss.
        for (e in 1:length(location[[i]])) {
            #test
            #e=1
            #test
            loss=var(outer_range[1:location[[i]][e]])*(location[[i]][e]-1)+var(outer_range[location[[i]][e]+1]:length(outer_range))*(length(outer_range)-location[[i]][e]-1)
            if (loss<best_loss) {
                best_loss=loss
                change_point[i+1]=location[[i]][e]
            }
            #cat("e",e,"\n")
            #cat("loss",loss,"\n")
            #cat("best_loss",best_loss,"\n")
            #cat("change_point",change_point,"\n")
        }
    }
    # Delete the first change point wich is fixed as 1
    change_point=change_point[-1]
    return(change_point)

}
RangeToPoint(x,n,location)
n=peakrange$n_peak_range
location=peakrange$peak_ranges

