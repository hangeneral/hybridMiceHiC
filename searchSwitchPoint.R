searchBelowLine <- function(inputVector, drawPlot=TRUE,...){
    inputVector <- sort(inputVector)
    inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
    slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
    xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum)
    #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
    y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
    if(drawPlot){  #if TRUE, draw the plot
        plot(1:length(inputVector), inputVector,type="l",...)
        b <- y_cutoff-(slope* xPt)
        abline(v= xPt,h= y_cutoff,lty=2,col=8)
        points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
        abline(coef=c(b,slope),col=2)
        title(paste("x=",xPt,", y=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x, Fold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
        #axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
    }
    return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
        yPt <- myVector[x]
        b <- yPt-(slope*x)
        xPts <- 1:length(myVector)
        # minimum number of points smaller than the line [x, yPt]
        return(sum(myVector<=(xPts*slope+b)))
}

searchAboveLine <- function(inputVector, drawPlot=TRUE,...){
    inputVector <- sort(inputVector)
    inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
    slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
    xPt <- floor(optimize(numPts_above_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum)
    #Find the x-axis point where a line passing through that point has the minimum number of points above it. (ie. tangent)
    y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
    if(drawPlot){  #if TRUE, draw the plot
        plot(1:length(inputVector), inputVector,type="l",...)
        b <- y_cutoff-(slope* xPt)
        abline(v= xPt,h= y_cutoff,lty=2,col=8)
        points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
        abline(coef=c(b,slope),col=2)
        title(paste("x=",xPt,", y=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x, Fold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
        #axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
    }
    return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#this is an accessory function, that determines the number of points above a diagnoal passing through [x,yPt]
numPts_above_line <- function(myVector,slope,x){
        yPt <- myVector[x]
        b <- yPt-(slope*x)
        xPts <- 1:length(myVector)
        # minimum number of points greater than the line [x, yPt]
        return(sum(myVector>=(xPts*slope+b)))
}




