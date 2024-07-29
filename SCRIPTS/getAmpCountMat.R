getAmpCountMat <- function(ampMat,amp_ids,norm_flag=NULL) {

	############## getAmpCountMat #################
	### Used to normalize by mean per Sample for a given amplicon-by-sample count matrix
	## norm_flag =1 normalize by mean of amplicon-RC across all the samples
	## norm_flag =2 normalize by mean of Sample-RC across all the amplicons
	## norm_flag =3  First,zero center by mean and then normalize by unit-vector-norm
	    message("   -- Normalizing by mean-coverage per sample")
        ampMat.norm = c()
        ind = rownames(ampMat) %in% amp_ids

        if (norm_flag==1) {
                ampMat.norm = apply(ampMat,1,function(x){x/mean(x)})
        }
        else if (norm_flag ==2) {
                ampMat.norm = apply(ampMat,2,function(x){x/mean(x)})
        }
        else if (norm_flag==3) {
                ampMat.std = apply(ampMat,2,function(x) {x -mean(x)})
                ampMat.std = apply(ampMat.std,2,function(x) {x/sqrt(sum(x^2))})
                ampMat.norm = ampMat.std
        }
        else {
                ampMat.norm = ampMat
        }
        geneAmpMat = ampMat.norm
        return(geneAmpMat)
  }
