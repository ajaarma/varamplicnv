getMixAmpRegFile <- function(bedFile,RegFile,skip_arg=2,gene_id) {
	
	ampData = read.table(bedFile,sep="\t",skip=skip_arg)
	ampData[,4] = gsub(" ","",ampData[,4])
	exon_list = rep(0,dim(ampData)[1])
	ampData = cbind(ampData,exon_list)

	
	regionData = read.table(regionFile,sep="\t",skip=skip_arg)
	regionData[,4] = gsub(" ","",regionData[,4])
	#geneIDs = unique(as.vector(regionData[,4]))
	geneIDs = gene_id
	
	count_exon = 0
	for (ele in geneIDs) {
		print(ele)
		#gene_id = ele
		ind = which(as.vector(regionData[,4])==ele)
		
		#for ( i in 1:dim(regionData)[1]) {
		chr_num = unique(as.vector(regionData[ind,1]))
		print(chr_num)
	
		ind_chr = which(as.vector(ampData[,1])==chr_num)
		tmp.amp.frame = ampData[ind_chr,]	
		
		count_exon = 0
		
		for ( i in ind) {
			print(i)
			reg_st = as.numeric(regionData[i,2])
			reg_en = as.numeric(regionData[i,3])
			diff_reg = reg_en - reg_st
			count_exon = count_exon+1
			print(count_exon)
		
			#for(j in 1:dim(ampData)[1]){
			for(j in 1:dim(tmp.amp.frame)[1]){
				amp_st = as.numeric(tmp.amp.frame[j,2])
				amp_en = as.numeric(tmp.amp.frame[j,3])
				diff_amp = amp_en - amp_st

				#cat(reg_st,"\t",reg_en,"\t",tmp.amp.frame[j,1],"\t",tmp.amp.frame[j,2],"\t",tmp.amp.frame[j,3],"\t",tmp.amp.frame[j,4],"\n",file="SKI_reg_amp.txt",append=TRUE)
		
				if ((amp_st <= reg_st) && (amp_en >= reg_st)){
					#cat(reg_st,"\t",reg_en,"\t",tmp.amp.frame[j,1],"\t",tmp.amp.frame[j,2],"\t",tmp.amp.frame[j,3],"\t",tmp.amp.frame[j,4],"\t",0,"\t",count_exon,"\t",j,"\n",file="SKI_reg_amp.txt",append=TRUE)
					
					tmp.amp.frame[j,7] = count_exon
				}
				if ((amp_st <=reg_en) && (amp_en >= reg_en)) {
					#cat(reg_st,"\t",reg_en,"\t",tmp.amp.frame[j,1],"\t",tmp.amp.frame[j,2],"\t",tmp.amp.frame[j,3],"\t",tmp.amp.frame[j,4],"\t",1,"\t",count_exon,"\t",j,"\n",file="SKI_reg_amp.txt",append=TRUE)
					tmp.amp.frame[j,7] = count_exon
				}
				
				if ((amp_st >=reg_st) && (amp_en <= reg_en)) {
					#cat(reg_st,"\t",reg_en,"\t",tmp.amp.frame[j,1],"\t",tmp.amp.frame[j,2],"\t",tmp.amp.frame[j,3],"\t",tmp.amp.frame[j,4],"\t",2,"\t",count_exon,"\t",j,"\n",file="SKI_reg_amp.txt",append=TRUE)
					tmp.amp.frame[j,7] = count_exon
				}

			}

		}

	}
	return(tmp.amp.frame)
}





getAmpID_BedReg <- function(bedFile,regionFile) {

	system("rm /mnt/test/NILS/BED/Reg_Amp_Intersect*")	
	#ampData = read.table(bedFile,sep="\t",skip=skip_arg)
	print("Here1")
	ampData = read.table(bedFile,sep="\t")
	
        ampData[,4] = gsub(" ","",ampData[,4])
        exon_list = rep(0,dim(ampData)[1])
        ampData = cbind(ampData,exon_list)


	print("After")
        #regionData = read.table(regionFile,sep="\t",skip=skip_arg)
        regionData = read.table(regionFile,sep="\t")
        regionData[,4] = gsub(" ","",regionData[,4])
        geneIDs = unique(as.vector(regionData[,4]))

	tmp.amp.reg.frame = list()

	count = 1
	for(ele in geneIDs) {
		
		ind = which(as.vector(regionData[,4]==ele))
		chr_num = unique(as.vector(regionData[ind,1]))
		
		ind_chr = which(as.vector(ampData[,1]==chr_num))
		tmp.amp.frame = ampData[ind_chr,]
		count_exon = 0
		
		print(ind)
		print(dim(tmp.amp.frame))
		for (i in ind) {
			print(i)
			reg_st = as.numeric(regionData[i,2])
			reg_en = as.numeric(regionData[i,3])
			diff_reg = reg_en - reg_st
			count_exon = count_exon+1	
			#for (j in 1:dim(ampData[ind_chr,])[1]) {
			for (j in 1:dim(tmp.amp.frame[ind_chr,])[1]) {
				amp_st = as.numeric(tmp.amp.frame[j,2])
				amp_en = as.numeric(tmp.amp.frame[j,3])
				diff_amp = amp_en - amp_st
				if ((amp_st <= reg_st) && (amp_en >= reg_st)){
                                        cat(chr_num,"\t",reg_st,"\t",reg_en,"\t",ele,"\t",as.vector(tmp.amp.frame[j,1]),"\t",tmp.amp.frame[j,2],"\t",tmp.amp.frame[j,3],"\t",tmp.amp.frame[j,4],"\t",0,"\t",count_exon,"\t",j,"\n",file="/mnt/test/NILS/BED/Reg_Amp_Intersect.txt",append=TRUE)
                                        tmp.amp.frame[j,7] = count_exon
                                }
                                if ((amp_st <=reg_en) && (amp_en >= reg_en)) {
                                        cat(chr_num,"\t",reg_st,"\t",reg_en,"\t",ele,"\t",as.vector(tmp.amp.frame[j,1]),"\t",tmp.amp.frame[j,2],"\t",tmp.amp.frame[j,3],"\t",tmp.amp.frame[j,4],"\t",1,"\t",count_exon,"\t",j,"\n",file="/mnt/test/NILS/BED/Reg_Amp_Intersect.txt",append=TRUE)
                                        tmp.amp.frame[j,7] = count_exon
                                }
                                
                                if ((amp_st >=reg_st) && (amp_en <= reg_en)) {
                                        cat(chr_num,"\t",reg_st,"\t",reg_en,"\t",ele,"\t",as.vector(tmp.amp.frame[j,1]),"\t",tmp.amp.frame[j,2],"\t",tmp.amp.frame[j,3],"\t",tmp.amp.frame[j,4],"\t",2,"\t",count_exon,"\t",j,"\n",file="/mnt/test/NILS/BED/Reg_Amp_Intersect.txt",append=TRUE)
                                        tmp.amp.frame[j,7] = count_exon
                                }

			}

		}
		tmp.amp.reg.frame[[count]] = tmp.amp.frame
		count = count+1

	}

	system(paste("sort -n -k1.4 -k6,6n /mnt/test/NILS/BED/Reg_Amp_Intersect.txt > /mnt/test/NILS/BED/tmp.unsort.txt; mv /mnt/test/NILS/BED/tmp.unsort.txt /mnt/test/NILS/BED/Reg_Amp_Intersect_sort.txt",sep=""))
	return(tmp.amp.reg.frame)
}
