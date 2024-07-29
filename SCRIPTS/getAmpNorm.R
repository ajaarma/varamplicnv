getAmpNorm <- function(ampMat,bedFile) {

	bedFile = read.csv(bedFile,sep="\t",skip=2,header=F)
	colnames(bedFile) = c("Chromosome","Start","End","Amplicon","BP","Strand")
	bedFile$diff = abs(bedFile[,2] - bedFile[,3])
	print(dim(ampMat))
	for ( i in 1:dim(ampMat)[1]) {
		
		amp_id = rownames(ampMat)[i]
		ind = which(bedFile$Amplicon==amp_id)
		val = bedFile$diff[ind]
		ampMat[i,] = ampMat[i,]/val
	}
	mat_list = list()
	mat_list[[1]] = bedFile
	mat_list[[2]] = ampMat
	return(mat_list)
}
