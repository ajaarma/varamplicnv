getGeneGCMatNorm <- function(ampMat,gcContent,amp_ids) {

        gcData <- read.csv(gcContent)
        tmp.mat = ampMat
        for ( i in 1:dim(tmp.mat)[1]) {
                ind = which(gcData$AMP_ID==rownames(tmp.mat)[i])
                gc_vals = gcData$GC_Count[ind]
                tmp.mat[i,] = tmp.mat[i,]/gc_vals
        }
        return(tmp.mat)
}


