getPCAAnalysis <- function(ampMat,j, k) {

        #ampMat = t(ampMat)
        #ampMat = apply(ampMat,2,function(x) x/mean(x))
        ampMat = rmZeroValRows(ampMat)
        ampMat_prep1 = t(apply(ampMat,1,function(x){x-mean(x)}))
        ampMat_prep = t(apply(ampMat_prep1,1,function(x){x/sqrt(sum(x^2))}))
        ampMat_new = ampMat_prep1

        cov_mat = cov(t(ampMat_new))
        #print("Before eigen")
        eigen_dat1 = eigen(cov_mat) #Eigen Value decompositon takes time hence it is once computed and saved as eigen_data.dat for whole matrix
        #save(eigen_dat,file="EIGEN/eigen_data3.dat") #eigen3.dat is computed in mean normalized data
        #save(eigen_dat,file="eigen_data4.dat") #eigen4.dat is computed in mean,GC corrected normalized data
        #save(eigen_dat,file="eigen_data5.dat") #eigen4.dat is computed in mean,GC corrected normalized data.No unit norm scaling
        #save(eigen_dat,file="eigen_data6.dat") #eigen6.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation.
        #save(eigen_dat,file="EIGEN/eigen_data7.dat") #GC & Mean normalized eigen7.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sample as observation.
        #print("After Eigen mat")
        #load("/mnt/test//mnt/test/SCRIPTS/cnv_panels/eigen_data.dat")
        print("After Saving Eigen Data")
        eigen_vec1 = eigen_dat1$vectors
        eigen_val1 = eigen_dat1$values
        print(eigen_val1[1:50])

        cum_sum = cumsum(eigen_val1)
        prop_var = cum_sum/sum(eigen_val1)

        print(prop_var[1:50])
        plot(prop_var[1:15],type="b",xlab="Number of First Principal Components",ylab="Proportion of Variance",main="Proportion of Variance Explained")
	x11()

        pca_ind = c(j:k)
	print(dim(eigen_vec1))
	print(dim(ampMat_new))
        x.hat.tmp = t(eigen_vec1[,pca_ind]) %*% ampMat_new
        #x.hat.tmp = t(eigen_vec[,pca_ind]) %*% ampMat
        x.hat = eigen_vec1[,pca_ind] %*% x.hat.tmp
        ampPCAMat = x.hat
        rownames(ampPCAMat) <- rownames(ampMat)
	return(ampPCAMat)
}


getPCAAnalysis2 <- function(ampMat, j,k) {
	#load("/mnt/test/SCRIPTS/cnv_panels/EIGEN/eigen_data137.dat")
	#load("/mnt/test/SCRIPTS/cnv_panels/EIGEN/eigen_data151.dat")
        #ampMat = t(ampMat)
        #ampMat = apply(ampMat,2,function(x) x/mean(x))
        ampMat = rmZeroValRows(ampMat)
        ampMat_prep1 = t(apply(ampMat,1,function(x){x-mean(x)}))
        ampMat_prep = t(apply(ampMat_prep1,1,function(x){x/sqrt(sum(x^2))}))
        ampMat_new = ampMat_prep1

        #cov_mat = cov(t(ampMat_new))
        print("Before eigen")
        #eigen_dat = eigen(cov_mat) #Eigen Value decompositon takes time hence it is once computed and saved as eigen_data.dat for whole matrix
        #save(eigen_dat,file="EIGEN/eigen_data3.dat") #eigen3.dat is computed in mean normalized data
        #save(eigen_dat,file="eigen_data4.dat") #eigen4.dat is computed in mean,GC corrected normalized data
        #save(eigen_dat,file="eigen_data5.dat") #eigen4.dat is computed in mean,GC corrected normalized data.No unit norm scaling
        #save(eigen_dat,file="eigen_data6.dat") #eigen6.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation.
        #save(eigen_dat,file="EIGEN/eigen_data9_mad.dat") #eigen9.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation.First 6 samples were removed.
        #save(eigen_dat,file="EIGEN/eigen_data11_mad.dat") #eigen9.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation.First 2 samples were removed.
        #save(eigen_dat,file="EIGEN/eigen_data12.dat") #eigen12.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Run137;
        #save(eigen_dat,file="EIGEN/eigen_data137.dat") #eigen12.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Run137;
        #save(eigen_dat,file="EIGEN/eigen_data137_mad.dat") #eigen12.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Run137;
        #save(eigen_dat,file="EIGEN/eigen_data151.dat") #eigen12.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Run151 after mad;
        #save(eigen_dat,file="EIGEN/eigen_data97.dat") #eigen97.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; 
        #save(eigen_dat,file="EIGEN/eigen_data97_rmSNP_Auto.dat") #eigen97.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Eigen after removing SNPs and duplicated amplicons; 5795 Amplicons;GC Normalization only Autosomes.
        #save(eigen_dat,file="EIGEN/eigen_data137_rmSNP_Auto.dat") #eigen97.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Eigen after removing SNPs and duplicated amplicons; 5795 Amplicons;GC Normalization only Autosomes.
        #save(eigen_dat,file="EIGEN/eigen_data151_rmSNP_Auto.dat") #eigen97.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Eigen after removing SNPs and duplicated amplicons; 5795 Amplicons;GC Normalization only Autosomes.
        #save(eigen_dat,file="EIGEN/eigen_data97_rmSNP_Auto_6013.dat") #eigen97.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Eigen after removing SNPs and duplicated amplicons; 5795 Amplicons;GC Normalization only Autosomes. 6013 Amplicons; No SNPs were removed.
        #save(eigen_dat,file="EIGEN/eigen_data82_rmSNP_Auto.dat") #eigen82.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Eigen after removing SNPs and duplicated amplicons; 5795 Amplicons;GC Normalization only Autosomes.
        #save(eigen_dat,file="EIGEN/eigen_data53_rmSNP_Auto.dat") #eigen82.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Eigen after removing SNPs and duplicated amplicons; 5795 Amplicons;GC Normalization only Autosomes.
        #save(eigen_dat,file="EIGEN/eigen_data97_Auto.dat") #eigen82.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Eigen after removing SNPs and duplicated amplicons; 5795 Amplicons;GC Normalization only Autosomes.
        #save(eigen_dat,file="EIGEN/eigen_data97_AutoDup.dat") #eigen82.dat is computed on mean centered and unit norm scaled. Each amplicon as random variable and sampke as observation. ALL Samples; Eigen after removing SNPs and duplicated amplicons; 5795 Amplicons;GC Normalization only Autosomes. Inclusion of redundant amplicons and their counts.
        print("After Eigen mat")
        #load("/mnt/test//mnt/test/SCRIPTS/cnv_panels/eigen_data.dat")
        print("After Saving Eigen Data")
        eigen_vec = eigen_dat$vectors
        eigen_val = eigen_dat$values
        print(eigen_val[1:50])

        cum_sum = cumsum(eigen_val)
        prop_var = cum_sum/sum(eigen_val)

        print(prop_var[1:50])
        plot(prop_var[1:100],type="b",xlab="Number of First Principal Components",ylab="Proportion of Variance",main="Proportion of Variance Explained")
	x11()

	print(dim(ampMat_new))
	print(dim(eigen_vec))
        pca_ind = c(j:k)
        x.hat.tmp = t(eigen_vec[,pca_ind]) %*% ampMat_new
        x.hat = eigen_vec[,pca_ind] %*% x.hat.tmp
        ampPCAMat = x.hat
        rownames(ampPCAMat) <- rownames(ampMat)
        return(ampPCAMat)
}



getFullPCAAnalysis2 <- function(ampMat.sorted,ampGCMat.sorted,j,k,p_flag = FALSE,auto_flag) {

        #ampPCAMat = getPCAAnalysis(ampGCMat.sorted,j,k)
	print(auto_flag)
	if(auto_flag=="AUTO") {
	        ampPCAMat = getPCAAnalysis2(ampGCMat.sorted,j,k)
	}
	if (auto_flag=="SEX") {
		ampPCAMat = getPCAAnalysis(ampGCMat.sorted,j,k)
	}
        #geneAmpPCAMat = ampPCAMat[ind_gene,]
	cat("Min :",min(ampPCAMat),"\n")
	cat("Max :",max(ampPCAMat),"\n")
        ampPCAMat.cov = addMeanCov(ampPCAMat,ampGCMat.sorted)
        #ampPCAMat.cov = ampPCAMat#,ampGCMat.sorted)
        #geneAmpPCAMat.cov = ampPCAMat.cov[ind_gene,]
        #geneAmpPCAMat.cov = ampPCAMat.cov
        #plotNormReadCountCustom(ampPCAMat.cov,3,j,k)
	if(p_flag==TRUE) {
        	plotNormReadCountCustom(ampPCAMat.cov,3,j,k)
        	#plotRatioCustom2(ampPCAMat.cov,1,j,k)
        	plotRatioCustom2(ampPCAMat.cov,1,j,k)
	}
	return(ampPCAMat.cov)
}

addMeanCov2 <- function(ampMat.sorted,ampGCMat.sorted){

	


}
