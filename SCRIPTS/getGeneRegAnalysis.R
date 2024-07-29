getGeneRegionAnalysis <- function(ind=c(1:100),ratio_flag=FALSE,j=3) {
	
	#ind = c(100:150)
	p_flag = FALSE
	ampPCAMat.cov = getFullPCAAnalysis2(ampMat.sorted,ampGCMat.sorted,j,k,p_flag)
	ampSVDMat.cov = getFullSVDAnalysis2(ampMat.sorted,ampGCMat.sorted,j,k,gcContent,p_flag)
	ampMDSMat.cov = getFullMDSAnalysis2(ampMat.sorted,ampGCMat.sorted,j,k,gcContent,p_flag)
	par(mfrow=c(2,2))
	if (ratio_flag==FALSE) {
		plotNormReadCountCustom(ampGCMat.sorted[ind_gene,][ind,],1,j,k)
		plotNormReadCountCustom(ampPCAMat.cov[ind_gene,][ind,],3,j,k)
		plotNormReadCountCustom(ampSVDMat.cov[ind_gene,][ind,],4,j,k)
		plotNormReadCountCustom(ampMDSMat.cov[ind_gene,][ind,],5,j,k)
	}
	else {

		plotRatioCustom2(ampGCMat.sorted[ind_gene,][ind,],3,j,k)
		plotRatioCustom2(ampPCAMat.cov[ind_gene,][ind,],1,j,k)
		plotRatioCustom2(ampSVDMat.cov[ind_gene,][ind,],5,j,k)
		plotRatioCustom2(ampMDSMat.cov[ind_gene,][ind,],6,j,k)
	}
}
