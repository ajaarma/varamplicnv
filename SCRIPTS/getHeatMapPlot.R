getHeatMapPlot <- function(geneMat,ind_gene,gene_id,my_pal_num=7,flag) {
	
	main_text = c()
	hclustfunc <- function(x) hclust(x, method="complete")
	distfunc <- function(x) dist(x,method="maximum")

	q_pal = round(my_pal_num/3)
	
	#colors = unique(c(seq(-3,-2,length=my_pal_num),seq(-2,0.5,length=my_pal_num),seq(0.5,2.5,length=my_pal_num)))
	#colors = unique(c(seq(min(geneMat),-0.5,length=my_pal_num),seq(-0.5,0.5,length=my_pal_num),seq(0.5,max(geneMat),length=my_pal_num)))
	##my_palette <- colorRampPalette(c("red", "white", "green"))(n = length(colors)-1)
	my_palette <- colorRampPalette(c("red", "white", "green"))(n = my_pal_num)
	
	if (flag == 1) {
                main_text = paste(" HeatMap Plot using PCA/FA for ",gene_id,sep="")
        }
        if (flag == 2) {
                main_text = paste("Plot of log2 Ratios for Standardization based Normalization method")
        }

        if (flag ==3) {
                main_text = c("HeatMap plot using GC Correct & Mean norm for ",gene_id,sep="")
        }

        if (flag==4) {
                main_text = paste("Plot of log2 Ratios for Mean and UNSC based Normalization method")
        }
        if (flag==5) {
                main_text <- paste("Heatmap Plot using SVD analysis for ",gene_id,sep="")
        }
        if (flag==6) {
                main_text <- paste("HeatMap plot using  SVD analysis for ",gene_id,sep="")
        }	

	#heatmap.2(t(geneMat[ind_gene,]),density.info="none",Colv=FALSE,trace="none",Rowv=TRUE,col=my_palette,main=main_text,margin=c(8,8),scale="col")
	heatmap.2(t(geneMat[ind_gene,]),density.info="none",Colv=FALSE,trace="none",Rowv=TRUE,col=my_palette,main=main_text,margin=c(8,8))
	#heatmap.2(t(geneMat[ind_gene,]),density.info="none",Colv=FALSE,trace="none",Rowv=TRUE,col=my_palette,main=main_text,margin=c(8,8),scale="row")
	#heatmap.2(t(geneMat[ind_gene,]),density.info="none",Colv=FALSE,trace="none",Rowv=TRUE,col=my_palette,main=main_text,margin=c(8,8),breaks=colors)

}
