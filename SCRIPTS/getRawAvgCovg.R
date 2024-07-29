getRawCovg <- function(ampMat,plot_path,repeat_flag=NULL) {

    message("-- Computing the mean of raw count and discarding the low covered samples < 100")	
	a1 = apply(ampMat,2,mean)
	#print(sort(a1))
	pdf(plot_path)
	#plot(a1,xlab="Samples",ylab="Raw average coverage",ylim=c(min(a1)-50,max(a1)+400),main="Raw average coverage per sample",pch=15,col="blue")
	plot(a1,xlab="Samples",ylab="Raw average coverage",ylim=c(50,max(a1)+400),main="Raw average coverage per sample",pch=15,col="blue",type="b")
	abline(h=100,col="red")
	text(c(1:length(a1)),a1,labels=names(a1),cex=0.7,adj=-0.1,srt=90)
	dev.off()

	### Filter out low quality samples having raw average coverage <100 ####
	a2 = a1[a1<100]
	message("   -- The Discarded low covered samples are: ",names(a2))
	ampMat = ampMat[,names(a1[a1>=100])]
	return(ampMat)	
}
