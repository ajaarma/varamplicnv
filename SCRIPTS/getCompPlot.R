getComparativePlot <- function(geneMatObj,range,flag=NULL,plot_flag=NULL,k=NULL){

        par(mfrow=c(1,length(geneMatObj)))
        if (flag==1){
                for (i in 1:length(geneMatObj)){
                        plotNormReadCountData(geneMatObj[[i]][range,],plot_flag[[i]],k)
                        #plotNormReadCountData(geneMatObj[[i]][range,],plot_flag[[i]])
                }
        }
        if (flag==2){
                for (i in 1:length(geneMatObj)){
                        plotRatioData(geneMatObj[[i]][range,],plot_flag[i],k)
                        #plotRatioData(geneMatObj[range,],plot_flag[[i]])
                }
        }
}


