plotNormReadCountData <- function(geneNormMat,flag=NULL,k=NULL) {

        print(dim(geneNormMat))

        #geneNormMat = log2(geneNormMat[75:150,])
        geneNormMat = geneNormMat#[75:125,]
        geneGCMat <- geneNormMat
        ### 490 36 dimension for FBN1 Gene #####
        main_text = c()
        if (flag==-1) {
                main_text = paste("Plot of Raw Count of Reads")
        }
        if (flag ==0) {
                        main_text = paste("Plot of Scaled Count of Reads normalized by Mean Amplicon")
                        #main_text = paste("Plot of Scaled Count of Reads normalized by Mean (Row/Column)")
                }
        if (flag==1) {
                        main_text = paste("Plot of Scaled count of Reads normalized by GC Content & Mean")
                }
        if (flag==2) {
                        main_text <- paste("Plot of Scaled count of reads")
                }
        if (flag==3) {

                        main_text <- paste("PCA/FA based normalization from 1 to",k," PC/Factor Loadings",sep="")
                }
                plot(c(1:dim(geneGCMat)[1]),geneGCMat[,1],type="l",ylim=c(min(geneGCMat),max(geneGCMat)),col=1,xlab="Amplicons",ylab="Unit Norm Count",main=paste(main_text))
                #plot(c(1:dim(geneGCMat)[1]),geneGCMat[,1],type="l",ylim=c(min(geneGCMat),max(geneGCMat)),col=1,xlab="Amplicons",ylab="Unit Norm Read Count",main=paste(main_text))
                #plot(c(1:dim(geneGCMat)[1]),geneGCMat[,1],type="l",ylim=c(-0.02,0.09),col=1,xlab="Amplicons",ylab="Normalized Read Count",main=paste(main_text))
                #plot(c(1:dim(geneGCMat)[1]),geneGCMat[,1],type="l",ylim=c(min(geneGCMat),max(geneGCMat)),col=1,xlab="Amplicons",ylab="Unit Norm Count",main=paste(main_text))

                for(i in 2:dim(geneGCMat)[2]) {
                        lines(c(1:dim(geneGCMat)[1]),geneGCMat[,i],col=i)
                }
}


