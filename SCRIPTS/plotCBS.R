suppressMessages(require("ggplot2"))
suppressMessages(require("gridExtra"))
options(bitmapType='cairo')

# gradient function
color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
      return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

#plotCBS <- function(segments,data,ext.data,amplicons,genes,sname,plot_dir) {
plotCBS <- function(segments,data,amplicons,genes,sname,plot_dir) {

    del_thresh <- -0.50
    dup_thresh <- 0.50
	# loop segments
	if (nrow(segments) == 0) {
		return;
	}
	out_list = list()

    print(head(segments))
	for (i in 1:nrow(segments)) {
		message(paste("  plotting segment ",i))
		# get segment position
		seg_start = segments$loc.start[i]
		seg_end = segments$loc.end[i]
		seg_chr = as.vector(segments$chrom[i])
		# select amplicons on relevant chr
		sub_chr <- data[data$chr == seg_chr,]
        # select exons overlapping CNV
        sub_gene <- genes[genes$chr == seg_chr & genes$start <= seg_end & genes$end >= seg_start, ]
        # Select all amplicons relevant to given gene and region
        gene_list = unique(sub_gene$gene)
        print(gene_list)

		#cnv_amp = sub_chr[sub_chr$start <= seg_end & sub_chr$end >= seg_start,]

        for(g  in 1:length(gene_list)) {

            gene_id = gene_list[g]
            seg_start = segments$loc.start[i]
            seg_end = segments$loc.end[i]

            cat(" -- Plotting for gene:",gene_id,"\n")
            #region_list = sub_gene[sub_gene$gene==gene_id,'region']
            region_list = genes[genes$gene==gene_id,'region']
            count = 1

            tmp_list = list()
            for (region_num in region_list) {
                tmp_amp = data[data$chr==seg_chr & data$gene==gene_id & data$region==region_num,]
                tmp_list[[count]] = tmp_amp
                count = count+1
            }
            region_amp = do.call("rbind",tmp_list)
            
            if (g ==1 & g < length(gene_list)) {
                print("First if")
                seg_start = seg_start
                tmp_end = max(region_amp$end)
                if (tmp_end <= seg_end) {
                    seg_end = tmp_end
                }
                else {
                    seg_end = seg_end
                }
            }
            else if (g >1 & g < length(gene_list)) {
                print("Second if")
                #seg_start = min(genes[genes$gene==gene_id,]$start)
                seg_start = min(region_amp$start)
                #seg_end = max(genes[genes$gene==gene_id,]$end)
                seg_end = max(region_amp$end)
            }
            else if (g == length(gene_list) & length(gene_list)>1) {
                print("Third if")
                #seg_start = min(gene[gene$gene==gene_id]$start)
                seg_start = min(region_amp$start)
                seg_end = seg_end
            }
            else{
                print("Fourth if")
                seg_start = seg_start
                seg_end = seg_end
            }
            
            cat(seg_start,"\t",seg_end,"\n")
            sub_gene <- genes[genes$gene == gene_id & genes$start <= seg_end & genes$end >= seg_start, ]
            region_list = sub_gene[sub_gene$gene==gene_id,'region']
            cnv_amp = sub_chr[sub_chr$start <= seg_end & sub_chr$end >= seg_start,]
            count = 1
            tmp_list = list()
            for (region_num in region_list) {
                tmp_amp = data[data$chr==seg_chr & data$gene==gene_id & data$region==region_num,]
                tmp_list[[count]] = tmp_amp
                count = count+1
            }
            sub_amp = do.call("rbind",tmp_list)
		    # select amplicons overlapping the CNV
		    #cnv_amp = sub_chr[sub_chr$start <= seg_end & sub_chr$end >= seg_start,]
 		    if (nrow(sub_amp) == 0) {
			    next	
		    }
		    # assign color codes (brown : out, blue : in and green : flanking)
            colnames(sub_amp)[2] = "sample"
		    sub_amp$group = 'ALL'
            sub_amp$color = 'brown'
            sub_amp[sub_amp$probe %in% cnv_amp$probe,'group']='FLANKING'
            sub_amp[sub_amp$probe %in% cnv_amp$probe,'color']='green'
		    sub_amp[sub_amp$start >= seg_start & sub_amp$end <= seg_end,'group'] = 'CNV'	
		    sub_amp[sub_amp$start >= seg_start & sub_amp$end <= seg_end,'color'] <- 'darkblue'
		    sub_amp[(sub_amp$group=="CNV"& sub_amp$sample < del_thresh/2),'color'] <- 'darkblue'#'orange'
		    sub_amp[(sub_amp$group=="CNV"& sub_amp$sample <= del_thresh),'color'] <- 'darkblue'
		    sub_amp[(sub_amp$group=="CNV" & sub_amp$sample > dup_thresh/2),'color'] <- 'darkblue' #'orange'
		    sub_amp[(sub_amp$group=="CNV" & sub_amp$sample >= dup_thresh),'color'] <- 'darkblue'
            #sub_amp$size = rep("30.005",dim(sub_amp)[1])

		    # scale :
		    min_x <- min(sub_amp$start)
		    max_x <- max(sub_amp$end)
		    plot_length <- max_x - min_x + 1
		    plot_width <- 1024
		    sub_amp$scaled_start <- plot_width*(sub_amp$start - min_x)/plot_length
		    sub_amp$scaled_end = plot_width*(sub_amp$end - min_x)/plot_length
		    sub_amp$scaled_centers = (sub_amp$scaled_start + sub_amp$scaled_end) / 2
		    # positions for box.
		    cnv_start_scaled = plot_width*(seg_start - min_x)/plot_length
		    cnv_end_scaled = plot_width*(seg_end - min_x)/plot_length
		    # limit y-scale.
		    sub_amp[which(sub_amp$sample < -2),'sample'] = jitter(-1.9,0.3)
		    sub_amp[which(sub_amp$sample > 2),'sample'] = jitter(1.9,0.3)
        
            message('Starting ggplot2')
		    g1 <- ggplot(sub_amp) + scale_y_continuous(limit=c(-2.9,2),labels=c('',-2,-1,0,1,2)) + scale_x_continuous(breaks=seq(0,1000,250),labels=c(format(min_x,big.mark=','), format(round(min_x + plot_length/(250/plot_width)),big.mark=','),format(round(min_x + plot_length/(500/plot_width)),big.mark=','),format(round(min_x + plot_length/(750/plot_width)),big.mark=','),format(round(min_x + plot_length/(1000/plot_width)),big.mark=',')))
		    # add the probe segments
		    #for (j in 1:nrow(sub_amp)) {
		    #    	g1 <- g1 + geom_rect(aes(xmin=scaled_start,ymin=(sample-0.01),xmax=scaled_end,ymax=(sample+0.01),fill=sub_amp[j,'color']),data=sub_amp[j,])
		    #}
		    g1 <- g1 + geom_rect(mapping=aes(xmin=scaled_start,ymin=(sample-0.01),xmax=scaled_end,ymax=(sample+0.01)),fill=sub_amp[,'color'],data=sub_amp,size=2.005)

		    # add boundaries
		    g1 <- g1 + annotate("rect",xmin=cnv_start_scaled,xmax=cnv_end_scaled,ymin=-2,ymax=2,color='red',alpha=.1)
		    # add titles
		    g1 <- g1 + labs(x='Genomic Position',y='Log2_Ratio',title=paste(sname,' : ',seg_chr,":",seg_start,"-",seg_end)) 
		    # add means
		    g1 <- g1 + geom_segment(aes(x=cnv_start_scaled,y=segments[i,'seg.mean'],xend=cnv_end_scaled,yend=segments[i,'seg.mean'],colour='CBS.Mean'),linetype = "dashed")
		    g1 <- g1 + geom_segment(aes(x=cnv_start_scaled,y=segments[i,'mean.nt.logR'],xend=cnv_end_scaled,yend=segments[i,'mean.nt.logR'],colour='Base.Mean'),linetype = "dashed")
		    
            # add gene
		    if (nrow(sub_gene) > 0) {
			    sub_gene$scaled_start <- plot_width*(sub_gene$start - min_x)/plot_length
			    sub_gene$scaled_end <- plot_width*(sub_gene$end - min_x)/plot_length
			    g1 <- g1 + geom_rect(aes(xmin=scaled_start,ymin=-2.3,xmax=scaled_end,ymax=-2.5),fill='black', data=sub_gene)
			    # print labels
			    g1 <- g1+annotate("text",x=(sub_gene$scaled_start+sub_gene$scaled_end)/2,y=-2.55,hjust=1,vjust=1,size=2.5,angle=45,label=sub_gene$region) 
			    # intron
			    g1 <- g1 + geom_segment(x=min(sub_gene$scaled_start),xend=max(sub_gene$scaled_end),y=-2.4,yend=-2.4,color='black')
			    # gene name
			    g1 <- g1 + annotate("text",x=(min(sub_gene$scaled_start)+max(sub_gene$scaled_end))/2,y=-2.15,label=sub_gene[1,'gene'])
		    }
		    
            message('saving')		
            #plot_dir = "/home/aakumar/CNV_PANEL/SCRIPTS/cnv_panels/OUTPUT/105/plots_single.undo_splits_by.0.5_single_1_9_21"
            #plot_dir = "/home/aakumar/CNV_PANEL/SCRIPTS/cnv_panels/OUTPUT/116/plots_single.undo_splits_by.0.5_single_5_9_24"
		    ggsave(paste(plot_dir,"/",paste(sname,"_",seg_chr,"_",seg_start,"_",seg_end,"_",gene_id,'.png',sep=''),sep=''),g1,width = 25,height=12, units='cm')
	        out_list[[i]] = list(sub_amp,sub_gene)
        }
        #print(i)
        message('plotting done')
    }
	return(out_list)
}

