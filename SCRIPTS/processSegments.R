mapSegCoord <- function(chr_num,seg_st_pos,seg_en_pos,amplidata) {

        ind_st = which(amplidata$chr==chr_num & amplidata$start==seg_st_pos)
        ind_en = which(amplidata$chr==chr_num & amplidata$start==seg_en_pos)
        ind_st_list = c(ind_st,ind_en)
        ind_low = min(ind_st_list)
        ind_max = max(ind_st_list)
        ind_range = c(ind_low:ind_max)
        pos_list = c()
        #for(ele in ind) {
        #       tmp_pos = amplidata$end[ele]
        #       pos_list = c(pos_list,tmp_pos)
        #}
        pos_list = amplidata$end[ind_range]
        end.pos.max = max(pos_list)
	end.pos.max.ind = which(amplidata$chr==chr_num & amplidata$end==end.pos.max)
	gene_id_list = unique(c(amplidata$gene[ind_st],amplidata$gene[end.pos.max.ind]))
	exon_list = unique(c(amplidata$region[ind_st],amplidata$region[end.pos.max.ind]))
	val_list = list()
	val_list[[1]] = end.pos.max
	val_list[[2]] = gene_id_list
	val_list[[3]] = exon_list
        #return(end.pos.max)
        return(val_list)
}

getSegmentsMap <- function(segments,data){

        pos_list = c()
        for (i in 1:dim(segments)[1]) {
                chr_num = segments$chrom[i]
                seg_st_pos = segments$loc.start[i]
                seg_en_pos = segments$loc.end[i]
                #ind_st = which(segments$loc.start==start_pos)
                #ind_en = which(segments$loc.start==end_pos)

                #end_pos_max = mapSegCoord(chr_num,end_pos,data)
                #end_pos_max = mapSegCoord(chr_num,seg_st_pos,seg_en_pos,data)
                val_list = mapSegCoord(chr_num,seg_st_pos,seg_en_pos,data)
		end_pos_max = val_list[[1]]
		gene_id_list = val_list[[2]]
		exon_list = val_list[[3]]
                #pos_list = c(pos_list,end_pos_max)
                segments$loc.end[i] <- end_pos_max
		segments$gene[i] <- paste(gene_id_list,sep=" ",collapse=";")
		segments$region[i] <- paste(exon_list,sep=" ",collapse="--")
		
        }
        #segments = cbind(segments,pos_list)
        return(segments)
}



processSegments <- function(segments.all,bed_file) {

		gene = vector(mode="character",length=dim(segments)[1])
                region = vector(mode="character",length=dim(segments)[1])
		segments = cbind(segments,gene,region)
		segments = getSegmentsMap(segments.all,bed_file)	

}
