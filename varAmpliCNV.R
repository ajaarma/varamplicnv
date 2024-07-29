### Creating Command-Line argument list
suppressMessages(require("optparse"))
suppressMessages(require("R.utils"))
sourceDirectory(paste("SCRIPTS/",sep=""))
#file.sources = list.files("SCRIPTS/",pattern="R$")
#sapply(paste("SCRIPTS/",file.sources,sep=""),source,.GlobalEnv)

option_list = list(
        make_option(c("-i","--ampCountMat"),type="character",
                    help="Amplicon count Matrix",metavar="character"),
        make_option(c("-o","--outDir"),type="character",
                    help="Output Directory path",metavar="character"),
        make_option(c("-b","--bedFile"),type="character",
                    help="Amplicon specific bedfile",metavar="character"),
        make_option(c("-r","--roi"),type="character",
                    help="Region Of Interest (ROI) file",metavar="character"),
        make_option(c("-c","--gcContent"),type="character",
                    help="GC Content",metavar="character"),
        make_option(c("-s","--gender"),type="character",
                    help="Gender of Input samples",metavar="character"),
        make_option(c("-a","--batch"),type="character",
                    help="Analysis Batch Name/Number",metavar="character"),
        make_option(c("-p","--propVarFrac"),default=0.80,type="numeric",
                    help="Proportion of variance to be removed. Default=80%",
                    metavar="numeric"),
        make_option(c("-n","--analysis"),default=0,type="numeric",
                    help="Direct Segmentation (DS=>0) or Amplicon 
                    Overlap Filtering (AOF=>1)",metavar="numeric")
    );
usage_arg = "
            Rscript varAmpliCNV.R 
                -i <sample-amplicon-count-Matrix> 
                -o <output-directory> 
                -b <amplicon-bed-file>
                -c <gc-content-file>
                -r <ROI-file>
                -s <gender-file>
                -a <batch Name/Number>
                -p <proprtion-of-variance for Auto/Sex: default:0.80>
                -n < DS=>0 or AOF=>1 Flag; default:0 >
            "

opt_parser = OptionParser(usage=usage_arg,option_list=option_list)
opt = parse_args(opt_parser)
#print(opt)
getCommandOptionCheck(opt,opt_parser)

#Extract input variables
ampCountData  = opt$ampCountMat
out_dir = opt$outDir
ampFile = opt$bedFile
regionFile = opt$roi
gcContent = opt$gcContent
genderFile = opt$gender
anal_type = opt$batch
pvFrac = opt$propVarFrac
ds_arg = opt$analysis

#Load packages
suppressMessages(require("Matrix"))
suppressMessages(require("graphics"))
suppressMessages(require("gplots"))
suppressMessages(require("DNAcopy"))
suppressMessages(require("gtools"))

#Create Output and Plots directories for each fo the Analysis Batches
dir_list = createFolders(out_dir,anal_type)
out_anal_dir <<- dir_list[[1]]
out_plot_dir <<- dir_list[[2]]

print(out_anal_dir)
print(out_plot_dir)

#system(paste("mkdir -p ",getwd(),"/OUTPUT/",anal_type,sep=""))
#system(paste("mkdir -p ",getwd(),"/PLOTS/",anal_type,sep=""))
###############################################################################
#
# Description: 1. Get Normalized, Unit norm scaled Amplicon-Sample Read Count 
#                matrix.
#	           2. The amplicons are sorted according to their coordinates	
#
#	Options: 1 : Mean Normalization across all the amplicons for a given sample
#		     2 : Mean Normalization and scaled to unit norm
#		     3 : No normalization. Raw count is intact
#		     4 : Zero mean centering and scaling to unit norm 
#
# Output: Sorted Amplicon-Sample Normalized count matrix according to given 
#         coordinates of Amplicons for Gene
#
################################################################################



# INITIALIZATIONS
###################################################

#DS Cut-Off
del_cf <<- -0.50
dup_cf <<-  0.50

message("\nAnalyzing Samples for the input Batch: ",anal_type,"\n")
cat(" -- Loading the raw mean coverage matrix\n ")
load(ampCountData)
#if (repeat_flag==0) {
	ampMat = getRawCovg(ampMat,paste(out_plot_dir,"/Raw_Avg_Covg",anal_type,
                                     ".pdf",sep=""))
	if(anal_type==82 | anal_type=='BATCH_2') {
	#ampMat = ampMat[,-c(1,13)]
		ampMat = ampMat[,-c(13,24)]
	#ampMat = ampMat[,-c(13)]
	#ampMat = ampMat
	}else if(anal_type==53 | anal_type=='BATCH_1') {
		ampMat = ampMat[,-c(18)]
	#ampMat = ampMat
	}else if (anal_type==105){
		ampMat = ampMat[,-c(12,16,18)]
	}else if(anal_type==136 | anal_type=='BATCH_4'){
	#ampMat = ampMat#[,-c(19)]
		ampMat = ampMat[,-c(16,17)] #removing s18 due to multiple CNV se
        #gments also has coverage 101 and s19 has mean coverage 109
	}else{
		ampMat = ampMat
	}

 #For batch run 82 - TAAD, sample s2 and S14 were removed because of multiple 
    #CNVs in same sample/gene
#ampMat = ampMat[,-c(18)] #For batch run  53 - TAAD, sample s25 is removed.
#ampMat = ampMat[,names(a1[a1>=90])] #For Deafness Batch run 253
#ampMat = ampMat[,-c(4, 12,13, 27, 28,29,30,36,37,38, 39, 40, 43, 44, 46)] #For 
#Sample ID 1009, removed due to low coverage,140841_s11.bam=72.63646; 140843_s13.bam=99.25997
########## GC Content Processing #######################
	gcData <- read.csv(gcContent)
	gcData[,4] = gsub(" ","",gcData[,4])

	ind = grep("MT",gcData[,4])
	if(length(ind) !=0) {
		gcData = gcData[-grep("^MT",gcData[,4]),]
	}else{
		gcData = gcData
	}
	a2 = rownames(ampMat)
	a1 = gcData$AMP_ID
	a3 = a2[match(a1,a2)]
	a11 = which(is.na(a3)==TRUE)
	if(length(a11)==0) {
		ampMat.sorted = ampMat[a3,] 
	}else {
		ampMat.sorted = ampMat[!is.na(a3),]
	}

	#ampFile = bedFile
	autoSome_flag = TRUE
	ampMat.auto = getAutosomalMat(ampMat.sorted,ampFile,autoSome_flag)[[1]]
	ampMat.sex.raw = getAutosomalMat(ampMat.sorted,ampFile,autoSome_flag)[[2]]
	ampMat.sex = getSampleSexCorrected(ampMat.sex.raw,ampMat.auto,genderFile)
	numSamples = dim(ampMat.sorted)[2] # Choice of Number of Prinicpal Compo
                 #nents to be incorporated for PCA based analysis.
	geneGCMat_list_auto = getGCMatAnalysis(ampMat.auto,ampFile,gcContent,
                                           anal_type,"AUTO")
	geneGCMat_list_sex = getGCMatAnalysis(ampMat.sex,ampFile,gcContent,
                                          anal_type,"SEX")
	ampGCMat.sorted.auto = geneGCMat_list_auto[[7]]
	ampGCMat.sorted.sex = geneGCMat_list_sex[[7]]
    save.image(file=paste(out_anal_dir,"/tmp2.RData",sep=""))

    #stop()
#### Computation of PCA/MDS on GC corrected and Mean 
#                         normalized data         ######################
    #cat("Entered values of j and sex_j are: ",j,"\t",sex_j,"\t",k,"\n")
	###################### MDS analysis for Autosomal targets ##########
    diff = 0
	mat_list_auto = getPCAViaMDSAnalysis(ampMat.auto,ampGCMat.sorted.auto,
                                         pvFrac,plot_flag=FALSE,anal_type,
                                         out_anal_dir,out_plot_dir,"AUTO")
	optPVIndex_auto = mat_list_auto[[1]]
	ampPCAMat.sorted.cov.auto = mat_list_auto[[4]]
	diff.auto = mean(diff)

	############### MDS analysis for SEX chromosomal targets ##############
	
    diff = 0
	mat_list_sex = getPCAViaMDSAnalysis(ampMat.sex,ampGCMat.sorted.sex,
                                        pvFrac,plot_flag=FALSE,anal_type,
                                        out_anal_dir,out_plot_dir,"SEX")

	optPVIndex_sex = mat_list_sex[[1]]
	ampPCAMat.sorted.cov.sex = mat_list_sex[[4]]
	diff.sex = mean(diff)
	
    cat("Number of PCs removed for autosomal and sex chromosome respectively: ",
        optPVIndex_auto-1,"\t",optPVIndex_sex-1,"\t",numSamples,"\n")

########### Computing the Log-R ratios for Leave-one-out normalization ########
	ampPCAMat.sorted.cov.full = rbind(ampPCAMat.sorted.cov.auto, 
                                      ampPCAMat.sorted.cov.sex)
    ampPCAMat.sorted.cov.auto.ratio = getLORatio(ampPCAMat.sorted.cov.auto,1,
                                            optPVIndex_auto,numSamples)*diff.auto
	ampPCAMat.sorted.cov.sex.ratio  = getLORatio(ampPCAMat.sorted.cov.sex,1,
                                            optPVIndex_sex,numSamples)*diff.sex
    ampMat.cov.ratio.mat = rbind(ampPCAMat.sorted.cov.auto.ratio, 
                                 ampPCAMat.sorted.cov.sex.ratio)
	save(ampMat.cov.ratio.mat,file=paste(out_anal_dir,"/ampMat.cov.ratio.mat_",
                                         optPVIndex_auto,"_",optPVIndex_sex,"_",
                                         numSamples,sep=""))
    #stop()

########## Segmentation using DS-approach  ###################################
	if (ds_arg==0) {
		mat_list_PCA_auto = getDNACopy(ampPCAMat.sorted.cov.auto.ratio,ampFile)
		mat_list_PCA_sex = getDNACopy(ampPCAMat.sorted.cov.sex.ratio,ampFile)
		ind_arg = c(4:(dim(ampMat.sorted)[2]+3))
		segment.smoothed.CNA.object.auto.PCA = getCBSAnalysis(mat_list_PCA_auto,
                                                        ampMat.sorted,ind_arg)
		segment.smoothed.CNA.object.sex.PCA = getCBSAnalysis(mat_list_PCA_sex,
                                                        ampMat.sorted,ind_arg)
		segment.auto = segments.summary(segment.smoothed.CNA.object.auto.PCA)
		segment.sex = segments.summary(segment.smoothed.CNA.object.sex.PCA)
		segments.all = rbind(segment.auto,segment.sex)
		write.table(segments.all,sep="\t",file=paste(out_anal_dir,
                                                "/DS_Segment_output_ALL",
                                                optPVIndex_auto,"_",
                                                optPVIndex_sex,".txt",
                                                sep=""),append=F,row.names=F,
                                                col.names=T,quote=F)
		####### Filtering out the segments by Amplicons <10 and having segmentation 
        ####### standard deviation >1 ####
		### Keeping only those segments whose segment mean lies between[-0.5,+0.5]
		ind1 = which(segments.all$seg.mean <= del_cf)
		ind2 = which(segments.all$seg.mean >= dup_cf)
		ind_cnv = c(ind1,ind2)
		segments.all1 = segments.all[ind_cnv,]
		ind = which(segments.all1$seg.sd <=1)
		segments.all1 = segments.all1[ind,]
		ind = which(segments.all1$num.mark >=10)
		segments.all1 = segments.all1[ind,]
		
        #### Processing segments and annotating them for Genes and Exons #####
		write.table(segments.all1,sep="\t",file=paste(out_anal_dir,
                                                "/DS_Segment_output_",
                                                optPVIndex_auto,"_",
                                                optPVIndex_sex,".txt",
                                                sep=""),append=F,row.names=F,
                                                col.names=T,quote=F)
	}else{
########## Segmentation using DS-AOF alogirthm ###############################
		#processAOF(anal_type,j,sex_j,bed_file,regionFile, ampMat.cov.auto.ratio2,
        message(" -- Computing AOF ")
		processAOF(ampMat.cov.ratio.mat,ampFile,regionFile,anal_type,
                   optPVIndex_auto,optPVIndex_sex,numSamples,abs(del_cf),
                   out_anal_dir)
    }

