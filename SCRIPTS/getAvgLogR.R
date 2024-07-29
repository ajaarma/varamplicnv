GetAvgLogR <- function(segments,data) {
    message(' - computing Avg-LogR for the segments')
    for (i in 1:nrow(segments)) {
        # cnv details
        cnv_chr <- as.character(segments[i,'chrom'])
        cnv_start <- segments[i,'loc.start']
        cnv_end <- segments[i,'loc.end']
        scaled_end <- cnv_end - cnv_start + 1
        #print(segments[i,])
        tmp <- matrix(0,ncol=2,nrow=(cnv_end - cnv_start + 1))
        #print(dim(tmp))
        if(dim(tmp)[1] !=1) {
            sub_data <- data[which(data$chr == cnv_chr & data$end >= cnv_start & data$start <= cnv_end),]
            #sub_data <- data[which(data$chr == cnv_chr & data$start >= cnv_start & data$start <= cnv_end),]
            if (nrow(sub_data) == 0) {
                message('sub data not found. should not be possible...')
                return
            }
            for (j in 1:nrow(sub_data)) {
                amp_start <- sub_data[j,'start'] - cnv_start + 1
                amp_end <- sub_data[j,'end'] - cnv_start + 1
                #amp_logR <- sub_data[j,'sample']
                amp_logR <- sub_data[j,'logR']
                if (amp_start < 1) {
                    amp_start <- 1
                }
                if (amp_end > scaled_end) {
                    amp_end <- scaled_end
                }
                tmp[amp_start:amp_end,1] <- tmp[amp_start:amp_end,1] + amp_logR
                tmp[amp_start:amp_end,2] <- tmp[amp_start:amp_end,2] + 1
            }
            tmp <- tmp[which(tmp[,2] > 0),]
            mean <- sum(tmp[,1]) / sum(tmp[,2])
            segments[i,'mean.nt.logR'] <- mean
            }
        }
    return(segments)

}

mapSegCoord <- function(chr_num,seg_st_pos,seg_en_pos,amplidata) {

    ind_st = which(amplidata$chr==chr_num & amplidata$start==seg_st_pos)
    ind_en = which(amplidata$chr==chr_num & amplidata$start==seg_en_pos)
    ind_st_list = c(ind_st,ind_en)
    ind_low = min(ind_st_list)
    ind_max = max(ind_st_list)
    ind_range = c(ind_low:ind_max)
    pos_list = c()
    #for(ele in ind) {
    #   tmp_pos = amplidata$end[ele]
    #   pos_list = c(pos_list,tmp_pos)
    #}
    pos_list = amplidata$end[ind_range]
    end.pos.max = max(pos_list)
    return(end.pos.max)
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
        end_pos_max = mapSegCoord(chr_num,seg_st_pos,seg_en_pos,data)
        #pos_list = c(pos_list,end_pos_max)
        segments$loc.end[i] <- end_pos_max
    }
    #segments = cbind(segments,pos_list)
    return(segments)
}

readGenesRegionInfo <- function(region_file) {

    genes <- read.table(region_file,as.is=T,header=F)
    #colnames(genes) <- c('chr','start','end','gene.region','score','strand')
    colnames(genes) <- c('chr','start','end','gene.region')
    # keep track of global gene start/end
    gene_start = list()
    gene_end = list()
    # split gene|region to two cols
    for (r in 1:nrow(genes)) {
        spl = strsplit(genes$gene.region[r],"|",fixed=TRUE)[[1]]
        genes$gene[r] = spl[1]
        genes$region[r] = spl[2]
        if (!genes$gene[r] %in% names(gene_start) ){
            gene_start[[ genes$gene[r] ]] <- genes$start[r]
        }
        else if ( gene_start[[ genes$gene[r] ]] > genes$start[r]) {
            gene_start[[ genes$gene[r] ]] <- genes$start[r]
        }
        if (!genes$gene[r] %in% names(gene_end)) {
            gene_end[[ genes$gene[r] ]] <- genes$end[r]
        }
        else if (gene_end[[ genes$gene[r] ]] < genes$end[r]) {
            gene_end[[ genes$gene[r] ]] <- genes$end[r]
        }
    }
    return(genes)
}

combineAmpGene <- function(amplicons,genes) {
    
    # add gene info.
    print(head(amplicons))
    for (r in 1:nrow(amplicons)) {
	#r=22338 #Case1
	#r=22322 #Case2
	#r=22348 #Case3
	#r=22356 #Case 4
	#r=22367 #Case 5
        #r_idx = which(genes$chr == amplicons$chr[r] & (genes$start - 500) <= amplicons$start[r] & (genes$end +500) >= amplicons$start[r])
        r_idx4 = which(genes$chr == amplicons$chr[r] & amplicons$start[r] <= genes$start & amplicons$end[r] >= genes$end) #Case4
        r_idx123_1 = which(genes$chr == amplicons$chr[r] & amplicons$start[r] >= genes$start & amplicons$start[r] <= genes$end) #case123
        r_idx123_2 = which(genes$chr == amplicons$chr[r] & amplicons$end[r] >= genes$start & amplicons$end[r] <= genes$end) #case123
        r_idx5 = which(genes$chr == amplicons$chr[r] & amplicons$start[r] >= genes$start & amplicons$end[r] <= genes$end) #case5
	    r_idx = unique(c(r_idx123_1,r_idx123_2,r_idx4,r_idx5))
	    if (length(r_idx)>1){
            if (length(r_idx)==3) {
		    cat(r,"\t",length(r_idx),"\t",r_idx,"\n")
            }
		    #r_idx = r_idx[1]
		    amplicons$gene[r] = paste(unique(genes$gene[r_idx]),collapse="|")
		    amplicons$region[r] = paste(genes$region[r_idx],collapse="|")
	    }else if (length(r_idx)!=0){
		#cat(r,"\t",length(r_idx),"\t",r_idx,"\n")
		    amplicons$gene[r] = genes$gene[r_idx]
            amplicons$region[r] = genes$region[r_idx]
	    }
        else if (length(r_idx)==0) {
		    #cat(r,"\t",length(r_idx),"\t",r_idx,"\n")
            amplicons$gene[r] = "NA"
            amplicons$region[r] = "NA"
        }
    }

    return(amplicons)
}

