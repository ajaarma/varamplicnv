getSampleSexCorrected <- function(ampMat.sex.raw,ampMat.auto.raw,gender_file) {

	sex.vector=NULL
	genderDat = read.table(gender_file,sep="\t")
	genderDat[,1] = gsub(" ","",genderDat[,1])
	genderDat[,2] = gsub(" ","",genderDat[,2])
	#head(genderDat)
		
	for(i in 1:dim(ampMat.sex.raw)[2]) {
		sample_names = colnames(ampMat.sex.raw)[i]
		#print(sample_names)
		ind = which(paste(as.vector(genderDat[,1]),".bam",sep="")==sample_names) #sample names with .bam extension
		#ind = which(paste(as.vector(genderDat[,1]),sep="")==sample_names)
		#print(ind)
		
		vals = as.vector(genderDat[ind,2])
		#print(vals)
		if (vals=="F"){
			sex.vector[i] <-0.5
		}
		else if (vals=="M") {
			sex.vector[i] <-1
		}
		else {
			
			sex.vector[i] <- sum(ampMat.sex.raw[,i])/ sum(ampMat.sex.raw[1:dim(ampMat.sex.raw)[1],] )*dim(ampMat.sex.raw)[2]
		}
	}

	#ampMat.auto2 = rmZeroValRows(ampMat.auto.raw)
	#propX = sum(ampMat.sex.raw[1:dim(ampMat.sex.raw)[1],])/dim(ampMat.sex.raw)[1]/sum(ampMat.auto2[1:dim(ampMat.auto2)[1],])*dim(ampMat.auto2)[1]

	ampMat.sex.raw = ampMat.sex.raw*sex.vector
	#print(sex.vector)

	return(ampMat.sex.raw)
}
