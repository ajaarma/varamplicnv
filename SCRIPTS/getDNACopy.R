getDNACopy <- function(ampMat,regionFile) {

	#regDat = read.table(regionFile,sep="\t",skip=2)
	regDat = read.table(regionFile,sep="\t")
	regDat[,1] = gsub(" ","",regDat[,1])
	regDat[,4] = gsub(" ","",regDat[,4])
	a12 = c()

	amp_common = intersect(as.vector(regDat[,4]),rownames(ampMat))
	ind_reg = as.vector(regDat[,4]) %in% amp_common
	regDat = regDat[ind_reg,]

	mat_list = list()
	mat_list[[1]] = regDat
	mat_list[[2]] = ampMat
	regDatNew = cbind(regDat[,c(4,1,2)],as.data.frame(as.matrix(ampMat)))
	regDatNew[,3] = regDatNew[,3]
	colnames(regDatNew) <- c("Clone","Chromosome", "Position",paste(colnames(ampMat),sep=","))
	mat_list[[3]] = regDatNew

	return(mat_list)	

}

getDNACopy2 <- function(ampMat,regionFile) {

        regDat = read.table(regionFile,sep="\t",skip=2)
        #regDat = read.table(regionFile,sep="\t")
        regDat[,1] = gsub(" ","",regDat[,1])
        #regDat[,8] = gsub(" ","",regDat[,8])
        regDat[,4] = gsub(" ","",regDat[,4])
        a12 = c()
        for ( i in 1:dim(regDat)[1]) {
                strs = strsplit(as.vector(regDat[i,1]),split="chr")[[1]][2]
                a12 = c(a12,strs)
        }

        a12[a12=="X"] <- 23
        a12[a12=="Y"] <- 24

        regDat1 = cbind(regDat,as.numeric(a12))
        regDat2 = regDat1[order(regDat1[,7]),]
        #regDat2 = regDat1[order(regDat1[,12]),]
        #regDat  = regDat2[,c(1:6)]
        regDat  = regDat2


        amp_common = intersect(as.vector(regDat[,4]),rownames(ampMat))
        #amp_common = intersect(as.vector(regDat[,8]),rownames(ampMat))
        #ind_reg = as.vector(regDat[,8]) %in% amp_common
        ind_reg = as.vector(regDat[,4]) %in% amp_common
        regDat = regDat[ind_reg,]

        ind_list = c()
        for (ele in regDat[,4]) {
        #for (ele in regDat[,8]) {
                ind = which(rownames(ampMat)==ele)
                ind_list = c(ind_list,ind)
        }
        ampMat.new = ampMat[ind_list,]
        mat_list = list()
        mat_list[[1]] = regDat
        mat_list[[2]] = ampMat.new
        regDatNew = cbind(regDat[,c(4,1,2)],as.data.frame(as.matrix(ampMat.new)))
        #regDatNew = cbind(regDat[,c(8,5,6)],as.data.frame(as.matrix(ampMat.new)))
        regDatNew[,3] = regDatNew[,3]+1
        colnames(regDatNew) <- c("Clone","Chromosome", "Position",paste(colnames(ampMat.new),sep=","))
        mat_list[[3]] = regDatNew
        #mat_list[[4]] = scale(regDatNew)

        #return(regDatNew)
        return(mat_list)

}


getCBSAnalysis <- function(mat_list_arg,ampMat_arg,ind_arg,chr_num=NULL) {

	regDatNew = mat_list_arg[[3]]
	print(ind_arg)
	print(dim(regDatNew))
	CNA.object <- CNA(regDatNew[,ind_arg], regDatNew$Chromosome,regDatNew$Position, data.type="logratio",sampleid=c(paste(colnames(regDatNew)[ind_arg],sep=",")))

	smoothed.CNA.object <- smooth.CNA(CNA.object)
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1,min.width=2,p.method="perm",alpha=0.01) #Run alternatively for 97 batch. 28.11.2017
	
	return(segment.smoothed.CNA.object)


}

plotCBSMethod <- function(meth1, meth2,ampMat_arg,ind_arg,chr_num) {

	par(mfrow=c(2,1))
	seg1 = getCBSAnalysis(meth1,ampMat_arg,ind_arg,chr_num)
	seg2 = getCBSAnalysis(meth2, ampMat_arg,ind_arg,chr_num)
}

getDNACopyAuto <- function(ampMat,amp_reg) {

	regDat = read.table(amp_reg,sep="\t")
        regDat[,1] = gsub(" ","",regDat[,1])
        regDat[,8] = gsub(" ","",regDat[,8])

	a12 = c()
        for ( i in 1:dim(regDat)[1]) {
                strs = strsplit(as.vector(regDat[i,1]),split="chr")[[1]][2]
                a12 = c(a12,strs)
        }

        a12[a12=="X"] <- 23
        a12[a12=="Y"] <- 24
		
	amp_id_list = c()
	chrom_list = c()
	pos_list = c()
	mat_list = list()
	
	for( i in 1:dim(ampMat)[1]) {
		ind = which(regDat[,8]==rownames(ampMat)[i])[1]
		chrom_num = regDat[ind,1]
		print(ind)
		print(regDat[ind,8])
		pos_num = regDat[ind,6]+1
		
		amp_id_list = c(amp_id_list,regDat[ind,8])
		chrom_list = c(chrom_list,chrom_num)
		pos_list = c(pos_list,pos_num)
	}
	#regDat.new = as.data.frame(cbind(amp_id_list,chrom_list,pos_list))
	regDat.new = cbind(cbind(amp_id_list,chrom_list,pos_list),as.data.frame(as.matrix(ampMat)))
	colnames(regDat.new) = c("Clone","Chromosome", "Position",paste(colnames(ampMat),sep=","))

	regDat.new$Position = as.numeric(regDat.new$Position)
	
	#mat_list[[1]] = regDat.new
	return(regDat.new)
	
}
		




