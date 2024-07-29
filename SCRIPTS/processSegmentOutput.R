processSegmentOutput <- function(segments.all,bedOverlapFile,mat_list_PCA,j,anal_type) {


	system(paste(" rm  OUTPUT/",anal_type,"/Segment_",anal_type,"_output_",j,".txt",sep=""))
	ovl_data = read.table(bedOverlapFile,sep="\t",header=TRUE,row.names=NULL)
	ovl_data[,1] = gsub(" ","",ovl_data[,1])
	ovl_data[,2] = gsub(" ","",ovl_data[,2])
	ovl_data[,3] = gsub(" ","",ovl_data[,3])
	print(ovl_data[1:5,1:5])	
	amp_names = mat_list_PCA[[3]]$Clone
	
	
		
	for (i in 1:dim(segments.all)[1]) {
		mark_st = as.vector(segments.all[i,3])
		mark_en = as.vector(segments.all[i,4])
		chr_num = as.vector(segments.all[i,2])
		print(chr_num)
		tmp.ovl = ovl_data[ovl_data$CHR==chr_num,]
		amp_st = tmp.ovl$AMP_ID[grep(mark_st,tmp.ovl$AMP_ST)]
		amp_en = tmp.ovl$AMP_ID[grep(mark_en,tmp.ovl$AMP_ST)]
		print(i)
		print(amp_st)
		print(amp_en)
		ind_list = c()
		for (e in amp_st ){
			tmp.ind = which(amp_names==e)
			ind_list = c(ind_list,tmp.ind)
		}
		for (e in amp_en) {
			tmp.ind = which(amp_names==e)
			ind_list = c(ind_list,tmp.ind)
		}
		ind_list = sort(ind_list)
		#cat(l_ind,h_ind,sep="\t")
		cat(ind_list)
		sink(file=paste("OUTPUT/",anal_type,"/Segment_",anal_type,"_output_",j,".txt",sep=""),append=T)
		if(length(ind_list) !=0) {
			range_list = c(ind_list[1]:ind_list[length(ind_list)])	
			#amp_list = amp_names[l_ind:h_ind]
			amp_list = amp_names[range_list]
			for (ele in amp_list) {
				#cat(segments.all[i,1],segments.all[i,2],segments.all[i,3],segments.all[i,4],segments.all[i,5],segments.all[i,6],segments.all[i,7],segments.all[i,8],segments.all[i,9],ele,"\n",sep="\t",file=paste("OUTPUT/97/Segment_",anal_type,"_output_",j,".txt",sep=""),append=T,row.names=F,col.names=T,quote=F)
				cat(segments.all[i,1],as.vector(segments.all[i,2]),segments.all[i,3],segments.all[i,4],segments.all[i,5],segments.all[i,6],segments.all[i,7],segments.all[i,8],segments.all[i,9],ele,"\n",sep="\t")
				#cat("\n",file=paste("OUTPUT/97/Segment_",anal_type,"_output_",j,".txt",sep=""),append=T,row.names=F,col.names=T,quote=F)

			}	
			cat("\n")
		}
		else {
			cat(segments.all[i,1],as.vector(segments.all[i,2]),segments.all[i,3],segments.all[i,4],segments.all[i,5],segments.all[i,6],segments.all[i,7],segments.all[i,8],segments.all[i,9],"OFF TARGET","\n",sep="\t")
		}
		sink()
		cat("\n")
	}
	cat("\n\n")
	
	
}
