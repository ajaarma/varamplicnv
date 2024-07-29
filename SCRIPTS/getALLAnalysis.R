getALLAnalysis <- function() {

        par(mfrow=c(3,2))
	p_flag=TRUE
	ampPCAMat.cov = getFullPCAAnalysis2(ampMat.sorted,ampGCMat.sorted,1,k,p_flag)
        ampPCAMat.cov = getFullPCAAnalysis2(ampMat.sorted,ampGCMat.sorted,2,k,p_flag)
        ampPCAMat.cov = getFullPCAAnalysis2(ampMat.sorted,ampGCMat.sorted,3,k,p_flag)

        x11()
#getFullPCAAnalysis2(ampMat.sorted,ampGCMat.sorted,k)

        par(mfrow=c(3,2))
	ampSVDMat.cov = getFullSVDAnalysis2(ampMat.sorted,ampGCMat.sorted,1,k,gcContent,p_flag)
        ampSVDMat.cov = getFullSVDAnalysis2(ampMat.sorted,ampGCMat.sorted,2,k,gcContent,p_flag)
        ampSVDMat.cov = getFullSVDAnalysis2(ampMat.sorted,ampGCMat.sorted,3,k,gcContent,p_flag)

        x11()
        par(mfrow=c(3,2))
	ampMDSMat.cov = getFullMDSAnalysis2(ampMat.sorted,ampGCMat.sorted,1,k,gcContent,p_flag)
	ampMDSMat.cov = getFullMDSAnalysis2(ampMat.sorted,ampGCMat.sorted,2,k,gcContent,p_flag)
	ampMDSMat.cov = getFullMDSAnalysis2(ampMat.sorted,ampGCMat.sorted,3,k,gcContent,p_flag)
}
