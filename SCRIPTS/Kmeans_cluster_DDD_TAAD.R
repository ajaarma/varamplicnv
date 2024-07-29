library("tidyverse")
library("cluster")
library("factoextra")
library("Ckmeans.1d.dp")
library("RColorBrewer")
k=3
result = Ckmeans.1d.dp(df$seg.mean,k)
colors = brewer.pal(k,"Dark2")
plot(df$seg.mean, col=colors[result$cluster], pch=result$cluster, cex=1.5,
     main="Optimal univariate clustering given k",
     sub=paste("Number of clusters given:", k))
abline(h=result$centers, col=colors, lty="dashed", lwd=2)
legend("bottomright", paste("Cluster", 1:k), col=colors, pch=1:k, cex=1, bty="n")

wss <- function(k){
  kmeans(df,k,nstart=10)$tot.withinss
}

#Segmentation DS - Blinfoded

taad_data = read.table("data/vcnv/figure_1_data.txt",sep="\t",header=T)
#taad_data = read.table("data/vcnv/batch_1/output/ds_segment_10amp_3_3.txt",sep="\t",header=T)

#Estimate number of clusters using Direct segmentation approach
ds_bic_list = c()
for (i in 1:dim(taad_data)[1]) {
  resTaad = Ckmeans.1d.dp(taad_data$seg.mean,k=i)
  ds_bic_list = c(ds_bic_list,resTaad$BIC)
}
k=which(ds_bic_list==max(ds_bic_list))
plot(ds_bic_list,type="b",col="blue",main="BIC values per choice of 'k' clusters",
     xlab="k clusters",ylab="BIC",
     ylim=c(min(ds_bic_list)-5,0)
     )

# Cluster plot for DS segmentation k=2

resTaad = Ckmeans.1d.dp(taad_data$seg.mean,k,y=taad_data$seg.sd)

plot(resTaad,xaxt="n",yaxt="n",lwd="3",xlab="DS-Segmentation Mean",ylab="Segment SD",
     main="Optimal weighted univariate k-means (k=2) clustering of TAAD data",
     cex=1.5,ylim=c(0,1.05),pch="16")
axis(1,at=seq(-1,1,by=0.1))
axis(2,at=seq(0,1.5,by=0.1))
abline(v=resTaad$centers,col=1:k,lty="dotted",lwd=2)
text(taad_data$seg.mean,taad_data$seg.sd,labels=taad_data$DID,cex=1.0,pos=3,font=1)

#Cluster plot for DS-AOF segmentation
aof_bic_list = c()
for (i in 1:dim(taad_data)[1]) {
  resAOF = Ckmeans.1d.dp(taad_data$mean.nt.log,k=i)
  aof_bic_list = c(aof_bic_list,resAOF$BIC)
}
plot(aof_bic_list,type="b",col="blue",main="BIC values per choice of 'k' clusters",
     xlab="k clusters",ylab="BIC",
     ylim=c(min(ds_bic_list)-30,0)
)
k=which(aof_bic_list==max(aof_bic_list))
resAOF = Ckmeans.1d.dp(taad_data$mean.nt.log,k,y=taad_data$seg.sd)
plot(resAOF,xaxt="n",yaxt="n",lwd="4",xlab="DS-AOF Segmentation Mean",ylab="Segment SD",
     main="Optimal weighted univariate k-means (k=3) clustering of TAAD data (DS-AOF)",
     ylim=c(0,1.05),pch="16")
axis(1,at=seq(-1,1,by=0.1))
axis(2,at=seq(0,1.5,by=0.1))
abline(v=resAOF$centers,col=1:k,lty="dotted",lwd=2)
#text(taad_data$mean.nt.log,taad_data$seg.sd,labels=taad_data$DID,cex=1.0,pos=3,font=1)

########## BT threshold ###############
st_data = read.table("data/vcnv/AllVarListPruned_rmBadSample_PC3_CutOff01_V4.3.csv",head=F,sep=",")
st_bic_list = c()
for (i in 1:dim(st_data)[1]) {
  resDST = Ckmeans.1d.dp(st_data$V10,k=i)
  st_bic_list = c(st_bic_list,resDST$BIC)
}
plot(st_bic_list,type="b",col="blue",main="BIC values per choice of 'k' clusters",
     xlab="k clusters",ylab="BIC",
     ylim=c(min(st_bic_list)-30,0)
)
k=which(st_bic_list==max(st_bic_list))
resDST = Ckmeans.1d.dp(st_data$V10,k,y=st_data$V7)
plot(resDST,xaxt="n",yaxt="n",lwd="4",xlab="ST Segmentation Mean",ylab="Segment SD",
     main="Optimal weighted univariate k-means (k=4) clustering of TAAD data (ST)",
     ylim=c(0,4.2),pch="16")
axis(1,at=seq(-1,1,by=0.1))
axis(2,at=seq(0,10,by=0.1))
abline(v=resDST$centers,col=1:k,lty="dotted",lwd=2)

#DOOF dataset
doof_data19 = read.table("data/vcnv/figure_2_data.txt",head=T,sep="\t")
doof_data = doof_data19[-19,]
doof_bic_list = c()
for (i in 1:dim(doof_data)[1]) {
  resDOOF = Ckmeans.1d.dp(doof_data$mean.nt.logR,k=i)
  doof_bic_list = c(doof_bic_list,resDOOF$BIC)
}
k=which(doof_bic_list==max(doof_bic_list))
plot(doof_bic_list,type="b",col="blue",main="BIC values per choice of 'k' clusters",
     xlab="k clusters",ylab="BIC",
     ylim=c(min(doof_bic_list)-10,0)
)
resDOOF = Ckmeans.1d.dp(doof_data$mean.nt.logR,k,y=doof_data$seg.sd)
plot(resDOOF,xaxt="n",yaxt="n",lwd="4",xlab="ST Segmentation Mean",ylab="Segment SD",
     main="Optimal weighted univariate k-means (k=2) clustering of DEAFNESS data (DS-AOF)",
     ylim=c(0,1.0),pch="16")
axis(1,at=seq(-1,1,by=0.1))
axis(2,at=seq(0,1.0,by=0.1))
abline(v=resDOOF$centers,col=1:k,lty="dotted",lwd=2)

