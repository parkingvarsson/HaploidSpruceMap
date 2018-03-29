#GT_data is the genotype data extracted from the filtered vcf file 
pca_input<-GT_data #GT_data is the genotype data extracted from the filtered vcf file (all samples)

new_pca<-rep(NA,dim(pca_input)[2]*dim(pca_input)[1])
dim(new_pca)<-dim(pca_input)


for (i in 1:dim(pca_input)[2]){
  new_pca[pca_input[,i]=="b",i]<-1
  new_pca[pca_input[,i]=="a",i]<-0
}
new_pca2<-as.data.frame(t(new_pca))

apply(new_pca2,2,mean,na.rm=T)->p.mean
for(i in 1:dim(new_pca2)[2]){
  new_pca2[is.na(new_pca2[,i]),i]<-p.mean[i]
}

prcomp(~.,data=new_pca2,na.action="na.exclude")->sample_pca
##summary(sample_pca)
plot(sample_pca$x[,1],jitter(sample_pca$x[,2]),pch=19,cex=0.5, xlab="PC 1 (10%)",ylab="PC 2 (7%)",las=1)
points(sample_pca$x[sample_pca$x[,2] > 5,1],jitter(sample_pca$x[sample_pca$x[,2] > 5,2]),col="red",pch=19,cex=0.5)
points(sample_pca$x[sample_pca$x[,2] < -5 & sample_pca$x[,1] > 0,1],jitter(sample_pca$x[sample_pca$x[,2] < -5 & sample_pca$x[,1] > 0,2]),col="blue",pch=19,cex=0.5)
points(sample_pca$x[sample_pca$x[,1] < -2,1],jitter(sample_pca$x[sample_pca$x[,1] < -2,2]),col="green",pch=19,cex=0.5)



## which individuals belong to the different cluster
c(1:1559)[sample_pca$x[,2] > 5]->cluster_1
c(1:1559)[sample_pca$x[,1] > 0 & sample_pca$x[,2] < -5]->cluster_2
c(1:1559)[sample_pca$x[,1] < -2]->cluster_3

