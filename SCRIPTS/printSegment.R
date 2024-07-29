printSegment <- function(segments.all,bedFile,amp_ids,anal_type,j,sex_j) {

	system(paste("rm OUTPUT/",anal_type,"/Segment_Amplicon*",sep=""))
	bedDat = read.table(bedFile,sep="\t")
	bedDat[,1] = gsub(" ","",bedDat[,1])
	bedDat[,2] = gsub(" ","",bedDat[,2])
	bedDat[,3] = gsub(" ","",bedDat[,3])
	bedDat[,4] = gsub(" ","",bedDat[,4])

	for (i  in 1:dim(segments.all)[1]) {
	
		print(segments.all[i,])	
		chr_num = as.vector(segments.all[i,2])
		st_coord1 = as.vector(segments.all[i,3])
		st_coord2 = as.vector(segments.all[i,4])
		tmp.frame = bedDat[bedDat[,1]==chr_num,]
		en_coord2 = as.vector(tmp.frame[which(tmp.frame[,2]==st_coord2),2])
		
		amp_1 = as.vector(bedDat[which(tmp.frame[,2]==st_coord1),4])
		amp_2 = as.vector(bedDat[which(tmp.frame[,2]==st_coord2),4])
	
		ind_1 = which(tmp.frame[,2]==st_coord1)
		if (length(ind_1) >1) {
			ind_1 = sort(ind_1,decreasing=T)
		}

		ind_2 = which(tmp.frame[,2]==st_coord2)
		if (length(ind_2) >1) {
			ind_2 = sort(ind_2,decreasing=T)

		}	

		bed_amp_list = as.vector(tmp.frame[,4])
		print(ind_1)
		print(ind_2)
		bed_amp_list = bed_amp_list[ind_1:ind_2]
		print(amp_1)
		print(amp_2)
		print(chr_num)
		#ind1 = which(amp_ids==amp_1)
		#ind2 = which(amp_ids==amp_2)
		
		#amp_id_list = amp_ids[ind1:ind2]
		#for (e in amp_id_list) {
		for (e in bed_amp_list) {
			cat(chr_num,"\t",st_coord1,"\t",en_coord2,"\t",e,"\n")
			#cat(chr_num,"\t",st_coord1,"\t",en_coord2,"\t",e,"\n",file=paste("OUTPUT/",anal_type,"/Segment_Amplicon_",j,"_",sex_j,".txt",sep=""),append=T,row.names=F,col.names=F,quote=F)
	
		}
		cat("\n",file=paste("OUTPUT/",anal_type,"/Segment_Amplicon_",j,"_",sex_j,".txt",sep=""),append=T,row.names=F,col.names=F,quote=F)
		
		bed_amp_list = c()	
		
	}

}
