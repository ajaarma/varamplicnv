getGeneAmpData <- function(bedFile,regionFile,st_coord, en_coord, chr_num,gene_id) {
#getGeneAmpData <- function(bedFile,gene_id,chr_num) {

        system(paste("python /mnt/test/SCRIPTS/cnv_panels/getCoordBed.py",bedFile,regionFile,chr_num,st_coord,en_coord,gene_id,sep=" "))
                #geneAmpData = read.table(paste("/mnt/test/SCRIPTS/cnv_panels/",chr_num,"_coord.txt",sep=""),sep="\t")
        geneAmpData = read.table(paste("/mnt/test/SCRIPTS/cnv_panels/",chr_num,"_",gene_id,"_coord.txt",sep=""),sep="\t")
		geneAmpData = geneAmpData[order(geneAmpData[,1]),]
        return(geneAmpData)
  }

getGeneAmpRegionData <- function(bedFile,regionFile,gene_id,chr_num,skip_arg=2) {
	
	ampData = read.table(bedFile,sep="\t",skip=skip_arg)	
	regionData = read.table(regionFile,sep="\t",skip=skip_arg)
	regionData[,4] = gsub(" ","",regionData[,4])
	gene.region.tmp = regionData[which(regionData[,4]==gene_id),]

	system(paste("python /mnt/test/SCRIPTS/cnv_panels/getCoordBed.py",bedFile,regionFile,chr_num,gene_id,sep=" "))
	geneAmpData = read.table(paste("/mnt/test/SCRIPTS/cnv_panels/",gene_id,"_coord.txt",sep=""),sep="\t")	
	return(geneAmpData)
}

getSortBedFile <- function(bedFile,skip_arg=2) {
	
	ampData = read.table(bedFile,sep="\t",skip=skip_arg)
	system(paste("sort -n -k1.4 -k2,2n ",bedFile, " > testBedFile.txt",sep=""))
}

getAmpDataCustom <- function(ampMat.raw,amp_reg_file) {

	ampRegDat  = read.table(amp_reg_file,sep="\t")
	ampIDs = unique(as.vector(ampRegDat[,8]))
	ampMat.raw = ampMat.raw[rownames(ampMat.raw) %in% ampIDs,]
	return(ampMat.raw)
}

getAutosomalMat <- function(ampMat.raw,amp_reg_file,autoSome_flag=TRUE) {

	ampRegDat = read.table(amp_reg_file,sep="\t")
	#ampRegDat = read.table(amp_reg_file,sep="\t",skip=2)
	ampRegDat[,4] = gsub(" ","",ampRegDat[,4])
	ampRegDat[,1] = gsub(" ","",ampRegDat[,1])
	mat_list = list()
	if(autoSome_flag) {
        message(" -- Processing Autosomal & Sex chromosome")
		ind1 = which(as.vector(ampRegDat[,1]) %in% c("chrX","X"))
		ind2 = which(as.vector(ampRegDat[,1]) %in% c("chrY","Y"))
		ind = c(ind1,ind2)
				
		tmp.reg.auto = ampRegDat[-ind,]
		tmp.reg.sex = ampRegDat[ind,]
		ampIDs.auto = unique(as.vector(tmp.reg.auto[,4]))
		ampIDs.sex = unique(as.vector(tmp.reg.sex[,4]))
		ampMat.raw.auto = ampMat.raw[rownames(ampMat.raw) %in% ampIDs.auto,]
		ampMat.raw.sex = ampMat.raw[rownames(ampMat.raw) %in% ampIDs.sex,]
	}
	else {
		ampIDs = unique(as.vector(ampRegDat[,4]))
		ampMat.raw = ampMat.raw[rownames(ampMat.raw) %in% ampIDs,]
	}

	mat_list[[1]] = ampMat.raw.auto
	mat_list[[2]] = ampMat.raw.sex
	
	return(mat_list)
}
