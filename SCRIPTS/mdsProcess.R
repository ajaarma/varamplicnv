rmZeroValRows <- function(ampMat_arg) {
	ampMat = ampMat_arg
    col_means = t(apply(ampMat,1,sum))
    tmp.mat = ampMat[which(col_means!=0),]
    return(tmp.mat)
}

mds <- function(mat,k=NULL) {
    mat.mean.row = t(apply(mat,1,function(x)x-mean(x)))
    mat.mean.col = apply(mat.mean.row,2,function(x)x-mean(x))
    mat.mean.zero = mat.mean.col
    X = mat.mean.zero
    return(X)
}

getMDS_XTXCustom <- function(X,pvFrac,out_anal_dir,plot_out_dir,anal_type=NULL,auto_sex=NULL) {
    
    message(' -- Computing MDS')

    #Initializations
    numSamples = dim(X)[2]
	X = rmZeroValRows(X)
	X = mds(X)
    M = t(X) %*% X
    eigen_xtx = eigen(M)
    vals = eigen_xtx$values
	eigen_vec = eigen_xtx$vectors
    eigen_val = eigen_xtx$values

    cum_sum = cumsum(eigen_val)
    prop_var = cum_sum/sum(eigen_val)

    #Computing optimal index based on proportion of variance (default = 0.80)
	if (auto_sex=="AUTO") {
		#ind2 = which.min(abs(prop_var -0.80))
		optPVIndex = which.min(abs(prop_var -pvFrac))
	}
	else {
		#ind2 = which.min(abs(prop_var -0.80))
		optPVIndex = which.min(abs(prop_var - pvFrac))
	}
  	if(optPVIndex==1) {
		optPVIndex = optPVIndex
	}
	else {
		optPVIndex = optPVIndex + 1
	}
	
    cat(prop_var,sep="\n",file=paste(out_anal_dir,"/Prop_var_",auto_sex,"_",
                                     optPVIndex,"_",numSamples,".txt",sep=""))
    pdf(paste(plot_out_dir,'/PC_component_',auto_sex,'_',numSamples,'.pdf',sep=''))
	plot(prop_var[1:numSamples],type="b",xlab="Number of First Principal Components",
         ylab="Proportion of Variance",main="Proportion of Variance Explained")
    
    dev.off()
	
    mds_ind = c(optPVIndex:numSamples)
    mix_mat = as.matrix(eigen_xtx$vectors[,mds_ind])
    mat.tmp1 = t(mix_mat) %*% t(X)
    mat.tmp = mix_mat %*% mat.tmp1
    rownames(mat.tmp) <- colnames(X)
	
    mat_list = list()
	mat_list[[1]] = optPVIndex
    mat_list[[2]] = numSamples
	mat_list[[3]] = t(mat.tmp)

    return(mat_list)
}

getPCAViaMDSAnalysis <- function(ampMat.sorted,ampGCMat.sorted,pvFrac,plot_flag = FALSE,
                                anal_type,out_anal_dir,plot_out_dir,auto_sex=NULL) {

    message(' -- Starting MDS analysis')

    mat_list = getMDS_XTXCustom(ampGCMat.sorted,pvFrac,out_anal_dir,plot_out_dir,anal_type,auto_sex)
	optPVIndex = mat_list[[1]]
    numSamples  = mat_list[[2]]
	ampGC_MDSMat.sorted = mat_list[[3]]
    
    #Add mean coverage back to MDS/PCA normalized matrix
    ampGC_MDSMat.sorted.cov = addMeanCov(ampGC_MDSMat.sorted,ampGCMat.sorted)
	
    if(plot_flag==TRUE) {
	    ##plotNormReadCountCustom(ampMDSMat.cov,5,j,k)
        #plotRatioCustom2(ampMDSMat.cov,6,j,k)
	    ##plotRatioCustom2(ampMDSMat.cov,6,j,k)
		mn = 1
	}

	mat_list[[4]] = ampGC_MDSMat.sorted.cov 
	#return(ampMDSMat.cov)
	return(mat_list)
}
