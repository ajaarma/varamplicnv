plotRatioData <- function(geneMat,flag=NULL,k=1) {

        main_text = c()

        if (flag == 1) {
                a11 = apply(geneMat,1,function(x){x/mean(x)})
                geneMat <- t(a11)
                main_text = paste("Plot of Ratios for Denoised Data using 1: ",k," PCA/FA",sep="")
        }
        if (flag == 2) {
                a11 = apply(geneMat,1,function(x) {x/mean(x)})
                geneMat <- t(a11)
                #geneMat  = geneMat
                main_text = paste("Plot of Ratios for Standardization based Normalization method")
                #main_text = c("Plot of Ratios for Mean normalization method")
                #main_text = c("Plot of Ratios for Mean and GC Normalized")
        }

        if (flag ==3) {
                a11 = apply(geneMat,1,function(x) {x/mean(x)})
                geneMat <- t(a11)
                #geneMat = log2(geneMat)
                #main_text = c("Plot of Log2 Ratios for Unormalized Raw read counts")
                #main_text = c("Plot of Log2 Ratios for Mean Normalized of Amplicons across Sample")
                main_text = c("Plot of Ratios for Mean and GC Normalized of Amplicons across Sample")
        }

        if (flag==4) {

                a11 = apply(geneMat,1,function(x) {x/mean(x)})
                a11[is.nan(a11)] <-0
                geneMat <- t(a11)
                #geneMat  = geneMat
                main_text = paste("Plot of Ratios for Mean and UNSC based Normalization method")

        }
        plot(c(1:dim(geneMat)[1]),geneMat[,1],type="l",ylim=c(min(geneMat),max(geneMat)),col=1,xlab="Amplicons",ylab="Ratio",main=paste(main_text))
        #plot(c(1:dim(geneMat)[1]),geneMat[,1],type="l",ylim=c(-5,2),col=1,xlab="Amplicons",ylab="Ratio",main=paste(main_text))
        #plot(c(1:dim(geneMat)[1]),geneMat[,1],type="l",ylim=c(-90,90),col=1,xlab="Amplicons",ylab="Ratio",main=paste(main_text))
        for(i in 2:dim(geneMat)[2]) {
                lines(c(1:dim(geneMat)[1]),geneMat[,i],col=i)
        }
}


