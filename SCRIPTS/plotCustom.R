plotRatioCustom2 <- function(geneMat,flag=NULL,j=NULL,k=36) {

    main_text = c()
	a11 = apply(geneMat,1,function(x)x/mean(x))
	a11[is.nan(a11)] <- 0
	geneMat = t(a11)
	geneMat = log2(geneMat)
	geneMat[is.nan(geneMat)] <- 0
	geneMat[is.infinite(geneMat)] <- 0
	
	print(min(geneMat))

    if (flag == 1) {
        main_text = paste("Plot of log2 Ratios for Denoised Data using ",j,":",k," PCA/FA on GC normalized data",sep="")
    }
    if (flag == 2) {
        main_text = paste("Plot of log2 Ratios for Standardization based Normalization method")
                #main_text = c("Plot of Ratios for Mean normalization method")
                #main_text = c("Plot of Ratios for Mean and GC Normalized")
    }
    if (flag ==3) {
        main_text = c("Plot of log2 Ratios for Mean and GC Normalized of Amplicons across Sample")
    }
    if (flag==4) {
        main_text = paste("Plot of log2 Ratios for Mean and UNSC based Normalization method")
    }
	if (flag==5) {
        main_text <- paste("Plot of log2 Ratios after SVD analysis ",j,":",k," components",sep="")
    }
	if (flag==6) {
                main_text <- paste("Plot of log2 Ratios after MDS analysis ",j,":",k," components",sep="")
    }
	
    plot(c(1:dim(geneMat)[1]),geneMat[,1],type="l",ylim=c(min(geneMat),max(geneMat)),col=1,xlab="Amplicons",ylab="Ratio",main=paste(main_text))
    
    for(i in 2:dim(geneMat)[2]) {
        lines(c(1:dim(geneMat)[1]),geneMat[,i],col=i)
    }

	return(geneMat)
}



plotNormReadCountCustom <- function(geneNormMat,flag=NULL,j=NULL,k=NULL) {

    geneNormMat = geneNormMat
    geneGCMat <- geneNormMat
    main_text = c()
    mode_text = "Mean"
    if (flag==-1) {
        main_text = paste("Plot of Raw Count of Reads")
    }
    if (flag ==0) {
        main_text = paste("Plot of Amplicons normalized by Mean ")
    }
    if (flag==1) {
        main_text = paste("Plot of Amplicons normalized by GC Content & Mean")
    }
    if (flag==2) {
        main_text <- paste("Plot of Amplicons of unormzalized count of reads")
    }
    if (flag==3) {
        main_text <- paste("PCA/FA based normalization from ", j," : ",k, " PC/Factor Loadings",sep="")
    }
	if (flag==4) {
		main_text <- paste("SVD analysis based normalization from ",j," : ",k," components",sep="")
	}
	if (flag==5) {
        main_text <- paste("MDS analysis based normalization from ",j," : ",k," components",sep="")
    }
	if (flag==6) {
        main_text <- paste("Log 2 Ratio plot of GC and Mean Normalized data after MAD criterion ",sep="")
    }
    plot(c(1:dim(geneGCMat)[1]),geneGCMat[,1],type="l",ylim=c(min(geneGCMat),max(geneGCMat)),col=1,xlab="Amplicons",ylab="Normalized Read Count",main=paste(main_text))
		
    for(i in 2:dim(geneGCMat)[2]) {
        lines(c(1:dim(geneGCMat)[1]),geneGCMat[,i],col=i)
    }
	
}

addMeanCov <- function(ampMat, ampGCMat) {

	if(min(ampMat) <0) {
		diff = abs(min(ampMat))+3
	}
	else {
		diff=0
	}
	cat("the Difference is: ",diff,"\n")
	ampGCMat = ampGCMat[rownames(ampMat),]
	mean_cov = apply(ampGCMat,2,mean)
	cat("The mean coverage to be added: \n")
    print(mean_cov)
	mean_cov = mean_cov+diff
	diff <<- mean_cov
	mean_rep = rep(mean_cov,dim(ampGCMat)[1])
	mean.mat.row = matrix(mean_rep,dim(ampGCMat)[1],dim(ampGCMat)[2],byrow=FALSE)
	ampMat = ampMat + mean.mat.row
	return(ampMat)
}

getDiff <- function(ampMat,ampGCMat) {

	if(min(ampMat) <0) {
        diff <- abs(min(ampMat))+3
    }
    else {
        diff<-0
    }
    cat("the Difference is: ",diff,"\n")
    ampGCMat = ampGCMat[rownames(ampMat),]
    mean_cov = apply(ampGCMat,2,mean)
    print(mean_cov)
    mean_cov = mean_cov+diff
	
    return(mean_cov)

}
addMeanCov2 <- function(ampMat, ampGCMat) {

	ampGCMat = ampGCMat[rownames(ampMat),]
	mean_cov = apply(ampGCMat,1,mean)
	mean_rep = rep(as.vector(mean_cov),dim(ampGCMat)[2])
	mean.mat.row = matrix(mean_rep,dim(ampGCMat)[1],dim(ampGCMat)[2],byrow=FALSE)
	ampMat = ampMat+mean.mat.row
	return(ampMat)

}


getLORatio <- function(ampMat,flag=1,j=1,k=36) {

	mat.tmp = mat.or.vec(dim(ampMat)[1],dim(ampMat)[2])
	rownames(mat.tmp) = rownames(ampMat)
	colnames(mat.tmp) = colnames(ampMat)

	
	for ( i in 1:dim(ampMat)[1]) {
		for (m in 1:dim(ampMat)[2]) {
			mat.tmp[i,m] = ampMat[i,m]/mean(ampMat[i,-m])
		}
	}

	ampMat.log = log2(mat.tmp)
	
    if (flag == 1) {
        main_text = paste("Plot of log2 Ratios for Denoised Data using ",j,":",k," PCA/FA on GC normalized data",sep="")
    }
    else if (flag == 2) {
        main_text = paste("Plot of log2 Ratios for Standardization based Normalization method")
    }

    else if (flag ==3) {
        main_text = c("Plot of log2 Ratios for Mean and GC Normalized of Amplicons across Sample")
    }
    else if (flag==4) {
        main_text = paste("Plot of log2 Ratios for Mean and UNSC based Normalization method")
	}
    else if (flag==5) {
        main_text <- paste("Plot of log2 Ratios after SVD analysis ",j,":",k," components",sep="")
    }
    else if (flag==6) {
        main_text <- paste("Plot of log2 Ratios after MDS analysis ",j,":",k," components",sep="")
    }
        
    #plot(c(1:dim(ampMat.log)[1]),ampMat.log[,1],type="l",ylim=c(-9,4),col=1,xlab="Amplicons",ylab="Ratio",main=paste(main_text))
	
    #for(i in 2:dim(ampMat.log)[2]) {
    #    lines(c(1:dim(ampMat.log)[1]),ampMat.log[,i],col=i)
    #}

    return(ampMat.log)
}

