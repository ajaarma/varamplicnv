getMDS_PCACustom <- function(X,k=NULL) {

	### Input matrix rows - amplicons and columns - Samples
	X = t(mds(X))
	print(dim(X))
        #T = dim(X)[2] -1
        #C = (1/T)*(X %*% t(X));
	C = cov(t(X))
        eigen_xxt = eigen(C)
	eigen_vec = eigen_xxt$vectors
        eigen_val = eigen_xxt$values
	eigen_val[eigen_val <0] <- 0
	save(eigen_xxt,file="EIGEN/eigen_xxt.dat")
        print(eigen_val[1:50])
        cum_sum = cumsum(eigen_val)
        prop_var = cum_sum/sum(eigen_val)

        print(prop_var[1:50])
        diag_mat = diag(1/sqrt(eigen_xxt$values))
	diag_mat[is.nan(diag_mat)] <-0
        whitenMat = diag_mat[,-k] %*% t(eigen_xxt$vectors[,-k])
        v2 = (1/sqrt(T))*whitenMat%*%X
        v3 = t(eigen_xxt$vectors)%*%X

	mds_pca_list = list()
	mds_pca_list[[1]] = X
	mds_pca_list[[2]] = v3
	
	return(mds_pca_list)
}
