processDupAmp <- function(ampMat,dupAmpFile) {


	a1 = read.table(dupAmpFile,sep="\t")
	a1_strs = strsplit(as.vector(a1[,2]),",")

	for (i in 1:length(a1_strs)) {

		amp_ids = a1_strs[[i]]
		tmp.mat = ampMat[rownames(ampMat) %in% amp_ids,]
		tmp.mat.sum = as.vector(apply(tmp.mat,1,sum))
		ind = which(tmp.mat.sum !=0)
		if (length(ind) !=0) {
			#cat(i,"\t",ind,"\n")
		#amp_id1 = rownames(tmp.mat)[ind]
		
			vals = ampMat[ind,]/length(amp_ids)
			ampMat[rownames(ampMat) %in% amp_ids,] <- vals
		}
	}
	return(ampMat)
}
