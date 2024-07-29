getPDFPlot <- function(mat,ind_gene,row_flag,method_text,pdf_file){ #Zero for row #One for Columns

        mat.tmp = mat
        pdf(file=pdf_file,width=22,height=18)
#ind1 = which(rownames(X) %in% rownames(geneAmpMat.PCA) ==TRUE)
        if(row_flag==0) {
                for(i in 1:dim(mat.tmp)[2]){
        #plot(mat.tmp[i,],main=paste("DOC plot with XTX based MDS for Sample ",i,sep=""),ylim=c(min(mat.tmp),max(mat.tmp)))
                        plot(mat.tmp[,i],main=paste("DOC plot with",method_text," for Sample ",i,sep=""),type="l",col=i,ylim=c(-1.2,0.9))
                        #points(ind1,mat.tmp[ind1,i],col="red")
                }
        }
        else {
                for(i in 1:dim(mat.tmp)[1]){
        #plot(mat.tmp[i,],main=paste("DOC plot with XTX based MDS for Sample ",i,sep=""),ylim=c(min(mat.tmp),max(mat.tmp)))
                        plot(mat.tmp[i,],main=paste("DOC plot with",method_text," for Sample ",i,sep=""))
                        points(ind1,mat.tmp[i,ind1],col="red")
                }
        }
        dev.off()
}
