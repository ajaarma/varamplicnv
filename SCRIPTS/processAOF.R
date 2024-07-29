require("DNAcopy")
require("ggplot2")
require("optparse")
require("gplots")
require("gtools")
options(bitmapType='cairo')
set.seed(42)
#file.sources = list.files("/home/aakumar/CNV_PANEL/SCRIPTS/V2/varamplicnv/SCRIPTS/",pattern="R$")
#sapply(paste("/home/aakumar/CNV_PANEL/SCRIPTS/V2/varamplicnv/SCRIPTS/",
#             file.sources,sep=""),source,.GlobalEnv)

processAOF <- function(ampMat.cov.ratio.mat,ampFile,regionFile,anal_type,
                       optPVIndex_auto,optPVIndex_sex,numSamples,thresh,
                       out_anal_dir){

    #out_path = "/home/aakumar/CNV_PANEL/SCRIPTS/cnv_panels/OUTPUT/136/"
    #amp_file = args[[3]] #"/home/aakumar/CNV_PANEL/DOOF/BED/testAmpRmSNPRmDupV8_sorted.bed"
    #ampFile = "/home/aakumar/CNV_PANEL/DOOF/BED/testAmpRmSNPRmDupV8_sorted.bed"
    #regionFile = "/home/aakumar/CNV_PANEL/DOOF/BED/EXONS/testSortGeneRegion_V8.bed"
    #i1 = "4_9"; j1="22"

    auto_sex_j = paste(optPVIndex_auto,optPVIndex_sex,sep='_')
    col_j = numSamples
    undo_type = anal_type

    del_thresh <- -1*thresh
    dup_thresh <-  1*thresh#0.50

    plot_dir = paste(out_anal_dir,"/plots_",anal_type,'_',thresh,"_",
                     auto_sex_j,"_",col_j,sep='')
    dir.create(plot_dir, showWarnings = FALSE)
    file.remove(file.path(plot_dir, list.files(plot_dir)))

    # gradient function
    color.gradient <- function(x, colors=c("red","yellow","green"), 
                               colsteps=100) {
      return( colorRampPalette(colors) (colsteps) [ 
             findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
    }

    #load(mat_path)
    
    if (exists("ampMat.cov.ratio.full")) {
        ampdata <- ampMat.cov.ratio.full
    } else if (exists("ampMat.cov.ratio.mat")){
        ampdata <- ampMat.cov.ratio.mat
    } else {
        message("needed matfile not found")
        ls()
        stop()
    }

    if (undo_type == '') {
        undo_type = 'none'
    }

    message("Reading design file")

    # read amplicon design file.
    amplicons <- read.table(ampFile,as.is=T,header=F,skip=0)
    colnames(amplicons) <- c('chr','start','end','name')
    rownames(amplicons) <- amplicons$name

    #intersect with amplicon-data matrix
    amplicons <- amplicons[rownames(ampdata),]

    #gene-Region-Of-Interest file
    message("Reading gene info");
    print(regionFile)
    genes <- readGenesRegionInfo(regionFile)

    #sort amplicons by chr/pos
    message("Sorting amplicons by Chromosomes/Positions");
    amplicons.sort <- amplicons[with(amplicons,order(amplicons$chr,amplicons$start)),]

    #combine Amplicons and Gene Information
    message("Combining Gene and Amplicon Information");
    amplicons.sort.genes <- combineAmpGene(amplicons.sort,genes)

    print('HERE')
    print(amplicons.sort.genes[1:5,])
    # construct SD-list of amplicons.
    amp.sds <- apply(ampdata,1,sd)

    out_list=list()

    for (s in 1:ncol(ampdata)){
        cat("\nThe iteration num: ",s,"\n")
        message(paste("Processing",colnames(ampdata)[s]))
        sample <- data.frame(sample=ampdata[,s],probe=rownames(ampdata),
                             stringsAsFactors=F)
        data <- merge(sample,amplicons.sort.genes, by.x=c("probe"),
                      by.y=c('name'),all=F)
        colnames(data) = c("probe","logR","chr","start","end","gene","region")
        
        message(paste(" - Segmenting sample ",s," (",colnames(ampdata)[s],") 
                      using AUTO/SEX PC components : ",optPVIndex_auto," and ",
                      optPVIndex_sex,sep=""))
        
        # set column types
        data$start <- as.numeric(data$start)
        data$logR <- as.numeric(data$logR)
        
        ## segment
        CNA.object <- CNA(data$logR,data$chr,data$start,data.type='logratio',
                          sampleid=colnames(ampdata)[s])
        smoothed.CNA.object <- smooth.CNA(CNA.object)
        segments.object <- segment(smoothed.CNA.object,verbose=1,min.width=2,
                                p.method='perm',alpha=0.01,undo.splits="none")
        segments.summ <- segments.summary(segments.object)
        segments.summ <- getSegmentsMap(segments.summ,amplicons.sort.genes)
        segments <- segments.summ

        message(' - filtering segments by number of amplicons >=10')
        segments.10 = segments[segments$num.mark >=10,]
        if (s==1) {
            write.table(segments.10,sep="\t",file=paste(out_anal_dir,"/",
                        "AOF_Segment_output_",thresh,"_",auto_sex_j,"_",
                        col_j,".Full.txt",sep=""),append=F,row.names=F,
                        col.names=T,quote=F)
        }else {
            write.table(segments.10,sep="\t",file=paste(out_anal_dir,"/",
                        "AOF_Segment_output_",thresh,"_",auto_sex_j,"_",
                        col_j,".Full.txt",sep=""),append=T,row.names=F,
                        col.names=F,quote=F)
        }

        message(' - filtering AOF segments by DS threshold DUP >= +0.38 
                and DEL <=-0.50')
        segments.10.DS = segments.10[(segments.10$seg.mean >=0.38)|
                                     (segments.10$seg.mean <= -0.28),]
        
        if(nrow(segments.10)>0) {
            segments.10.avg = GetAvgLogR(segments.10,data)
        }else {
            segments.10$mean.nt.logR <- numeric(0)
            segments.10.avg = segments.10
        }

        #message('-filtering segments based on threshold value')
        segments.10.avg.thresh = segments.10.avg[
                            (segments.10.avg$mean.nt.logR >=0.38) | 
                            (segments.10.avg$mean.nt.logR <=-0.20),]
        if (s ==1) {
            write.table(segments.10.avg.thresh,sep="\t",
                        file=paste(out_anal_dir,"/","AOF_Segment_output_",
                        thresh,"_",auto_sex_j,"_",col_j,".Filtered.txt",
                        sep=""), append=F,row.names=F,col.names=T,quote=F)
        } else {
            write.table(segments.10.avg.thresh,sep="\t",
                        file=paste(out_anal_dir,"/","AOF_Segment_output_",
                        thresh,"_",auto_sex_j,"_",col_j,".Filtered.txt",
                        sep=""),append=T,row.names=F,col.names=F,quote=F)
        }

        if (nrow(segments.10.avg.thresh)==0) {
            message(paste("No CNVs found in ",colnames(ampdata)[s]))
            next
        }
        out_list[[s]] = plotCBS(segments.10.avg.thresh,data,
                                amplicons.sort.genes,genes,
                                colnames(ampdata)[s],plot_dir)
    }

}
