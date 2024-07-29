getKernelPCACustom <- function(X,k=NULL,sqSigma=NULL){

        if(is.null(sqSigma)){
                sqSigma = 1
        }

        T = dim(X)[2]
        dist = t(X) %*% X
        sqNorms = as.matrix(diag(dist))
        dist = -2*dist
        a1 = rep(sqNorms,dim(dist)[2])
        a1.mat.row = matrix(a1,4842,4842,byrow=FALSE)
        a1.mat.col = matrix(a1,4842,4842,byrow=TRUE)

        dist1 = dist+a1.mat.row
        dist2 = dist1+a1.mat.col
        #dist = -exp(-1/(2*sqSigma)*dist)
        dist3 = -exp((-1/2*2)*dist2)

        eigen_pca = eigen(dist3)
        eigen_vec = eigen_pca$vectors
        a2 = length(eigen_pca$values)
        ind = sort(c(1:a2),decreasing=TRUE)
        eigen_vec_sort = eigen_vec[,ind]
        x.tmp = X %*% eigen_vec_sort[,1:36]
        x.hat = eigen_vec_sort[,1:36] %*% x.tmp
	return(x.hat)
}


