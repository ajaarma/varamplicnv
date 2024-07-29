extractGCMat <- function(ampMat,gcContent,amp_ids=NULL,anal_type=NULL,auto_sex=NULL) {

	    message(" -- Normalizing data for GC content:",auto_sex)
        gcData <- read.csv(gcContent)

        ampMat_rows = rownames(ampMat)
        ind_gc = gcData$AMP_ID %in% ampMat_rows
        ind = ampMat_rows %in% gcData$AMP_ID

        gc_count_frame = gcData[ind_gc,]
        ampGC_MAT = ampMat[ind,] #GC content matrix with same number of amplicons as amplicon Sample count matrix (ampMat)
        data_list = list() # Data list containing all the objects after GC correction
        data_list[[1]] = gc_count_frame # GC counts per amplicon coordinate
        data_list[[2]] = ampGC_MAT #read count matrix containing amplicon same as GC content data

        ## Currently the loess function works for R=3.2.1 and crashes for R=3.2.2
        GC_FRAC = gc_count_frame$GC_Count/(gc_count_frame$CHR_ST+1 - gc_count_frame$CHR_EN)
        loess_rel = loess(as.matrix(ampGC_MAT) ~ as.matrix(GC_FRAC))
        loess_predict = data.frame(predict(loess_rel,GC_FRAC))

        loess_factor = median(as.matrix(ampGC_MAT))/loess_predict
        names(loess_factor)[1] <- "Factor"
        ampGC_MAT.norm = round(as.matrix(ampGC_MAT)*loess_factor$Factor)

        data_list[[3]] <- ampGC_MAT.norm #Normalized amplicon read count after GC correction
        data_list[[4]] <- loess_rel # Loess regression output
        data_list[[5]] <- loess_predict #Loess regression output
        data_list[[6]] <- GC_FRAC #Loess regression output
	
        #Correlation plot for GC fraction and amplicon Read counts
        cor_list = c()
        for (i in 1:dim(ampGC_MAT)[2]) {
                cor_list = c(cor_list,cor(as.vector(GC_FRAC),as.vector(ampGC_MAT[,i])))
        }
	    #pdf(paste(getwd(),"/PLOTS/",anal_type,"/GC_Mat_plot_",auto_sex,".pdf",sep=""))
	    pdf(paste(out_plot_dir,"/GC_Mat_plot_",auto_sex,".pdf",sep=""))
        names(cor_list) <- colnames(ampGC_MAT)
        plot(cor_list,main="Correlation Plot between GC Content fraction and Read depth",xlab="Samples",ylab="Correlation Coefficient",ylim=c(-1,1),pch=16,type="b",col="blue")
        abline(h=max(cor_list),col="red")
        abline(h=min(cor_list),col="red")
	    dev.off()
	    #x11()
        return(data_list)
}

getGCMatAnalysis <- function(ampMat.sorted,bedFile,gcContent,anal_type=NULL,auto_sex=NULL) {

        ampIDs = unique(as.vector(rownames(ampMat.sorted)))
        geneGCMat_list <- extractGCMat(ampMat.sorted,gcContent,ampIDs,anal_type,auto_sex)#,geneAmpIDs) #Extract GC content matrix list
        geneGCMat <- getAmpCountMat(geneGCMat_list[[3]],ampIDs,2) #Mean normalized of amplicons read count data after GC correction
        geneGCMat_list[[7]] <- geneGCMat #GC content and mean normalized data 
        geneGCMat = rmZeroValRows(geneGCMat)
        return(geneGCMat_list)

}




