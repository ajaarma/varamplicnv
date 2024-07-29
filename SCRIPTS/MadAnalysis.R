getMADAnalysis <- function(ampMat.raw) {

	ampMat.raw.orig = ampMat.raw
	ampMat.raw.200 = rmLowCovAmp(ampMat.raw)
	ampMat.raw = ampMat.raw.200
	#ampMat.mean.norm = apply(ampMat.raw,2,function(x)x/mean(x))
	#a1.mad.sorted = sort(apply(ampMat.mean.norm,1,mad),decreasing=TRUE)
	#a1.mad.sorted = apply(ampMat.mean.norm,1,mad)
	a1.mad.sorted = apply(ampMat.raw,1,mad)
	print(sort(a1.mad.sorted,decreasing=TRUE)[1:50])
	a1.mad.sorted.02 = a1.mad.sorted[!a1.mad.sorted>=15000]
	names.a1.mad.sorted.02 = names(a1.mad.sorted.02)
	ampMat.raw = ampMat.raw.orig[names.a1.mad.sorted.02,]	
	return(ampMat.raw)	
	#return(a1.mad.sorted)	
}


getMADAnalysis2 <- function(ampMat.arg) {
	
	amp_list = rownames(ampMat.arg)
	ampMat.arg.orig = ampMat.arg
	#ampMat.arg = apply(ampMat.arg,2,function(x)x-mean(x))
	amp_union = c()	
	for (i in 1:dim(ampMat.arg)[2]) {
		x1 = ampMat.arg[,i] #- mean(ampMat.arg[,i])
		mad_score = mad(x1)
		upper_lim = median(x1)+3*mad_score
		lower_lim = median(x1)-3*mad_score
		x1.upper.amp = 	x1[x1>upper_lim]
		x1.lower.amp = x1[x1<lower_lim]
		amp_x = unique(c(names(x1.upper.amp),names(x1.lower.amp)))
		#amp_list = intersect(amp_list,amp_x)
		amp_list = union(amp_union,amp_x)
		#print(amp_list)
		#amp_list_new = intersect(amp_list)
	}
	ampMat.arg = ampMat.arg[!rownames(ampMat.arg) %in% amp_list,]
	return(ampMat.arg)
}


rmLowCovAmp <- function(ampMat.raw){
	
	a1.mean.cov = apply(ampMat.raw,1,mean)
	a1.mean.cov.200 = a1.mean.cov[a1.mean.cov<=200]
	a1.mean.cov.200.mat = ampMat.raw[names(a1.mean.cov.200),]
	a1.mean.cov.mat = ampMat.raw[!rownames(ampMat.raw) %in% names(a1.mean.cov.200),]
	#return(a1.mean.cov.200.mat)
	return(a1.mean.cov.mat)
}


rmBadSample <- function(ampMat.raw,ind) {
	
	mat = ampMat[,-ind]

}
getDensityPlot <- function(ampMat_arg,ind=c(1:10),main_text=NULL){
	
	ampMat.tmp = ampMat_arg[,ind]
	dens.10 = apply(ampMat.tmp,2,density)
	plot(NA, xlim=range(sapply(dens.10, "[", "x")), ylim=range(sapply(dens.10, "[", "y")),xlab="Read Depth",ylab="Density",main = main_text)
	mapply(lines, dens.10, col=1:length(dens.10))
	legend("topright", legend=names(dens.10), fill=1:length(dens.10))
}



