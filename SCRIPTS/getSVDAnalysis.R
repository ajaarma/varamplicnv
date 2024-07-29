getSVDAnalysis <- function(ampMat,j=NULL,k=NULL,gcContent=NULL,ampIDs = NULL) {
        ampMat = rmZeroValRows(ampMat)
        #ampMat_prep1 = apply(ampMat,2,function(x){x/mean(x)})
        ampMat.svd = svd(t(ampMat))
        svd_ind = c(j:k)
	print(svd_ind)
        #a1 = ampMat.svd$v[,-svd_ind] %*% diag(ampMat.svd$d)[-svd_ind,]
	print(ampMat.svd$d[1:50])
        a1 = ampMat.svd$v[,svd_ind] %*% diag(ampMat.svd$d)[svd_ind,]
        ampMat.norm = a1 %*% ampMat.svd$u
        rownames(ampMat.norm) <- rownames(ampMat)
        colnames(ampMat.norm) <- colnames(ampMat)
        #cov_mat = 
        ampMat.sort = c()
        ampMat.orig.sort = c()
        if(! is.null(gcData)){
                gcData = read.csv(gcContent)
                ind_list = c()
                for (i in 1:dim(gcData)[1]){
                        amp_id = gcData$AMP_ID[i]
                        ind = which(rownames(ampMat.norm)==amp_id)
                        if (length(ind)!=0) {
                                ind_list = c(ind_list,ind)
                        }
                }
                ampMat.sort = ampMat.norm[ind_list,]
                ampMat.orig.sort = ampMat[ind_list,]
        }
        ampMat_list = list()
        ampMat_list[[1]] = ampMat
        ampMat_list[[2]] = ampMat.norm
        ampMat_list[[3]] = ampMat.sort
        ampMat_list[[4]] = ampMat.orig.sort
        return(ampMat_list)
}

getFullSVDAnalysis2 <- function(ampMat.sorted,ampGCMat.sorted,j,k,gcContent,p_flag = FALSE) {

        ampGC_SVDMat_list = getSVDAnalysis(ampGCMat.sorted,j,k,gcContent)
        ampGC_SVDMat.sorted = ampGC_SVDMat_list[[3]]
        #geneAmpSVDMat = ampGC_SVDMat.sorted[ind_gene,]
        ampSVDMat.cov = addMeanCov(ampGC_SVDMat.sorted,ampGCMat.sorted)
        #geneAmpSVDMat.cov = ampSVDMat.cov[ind_gene,]
        #plotNormReadCountCustom(ampSVDMat.cov,4,j,k)	
	if (p_flag==TRUE) {
	        plotNormReadCountCustom(ampSVDMat.cov,4,j,k)
        #plotRatioCustom2(ampSVDMat.cov,5,j,k)
        	plotRatioCustom2(ampSVDMat.cov,5,j,k)
	}
	return(ampSVDMat.cov)

}

