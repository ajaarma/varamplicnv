sortMatByAmp <- function(ampMat,gcContent) {

        ind_list = c()
        gcData = read.csv(gcContent)
        for (i in 1:dim(gcData)[1]){
                amp_id = gcData$AMP_ID[i]
                ind = which(rownames(ampMat)==amp_id)
                if (length(ind)!=0) {
                        ind_list = c(ind_list,ind)
                }
        }
        ampMat.sort = ampMat[ind_list,]
        sort_obj = list()
        sort_obj[[1]] = ind_list
        sort_obj[[2]] = ampMat.sort

        return(sort_obj)

}
