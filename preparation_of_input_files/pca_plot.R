pca_input<-read.table("pca_analysis.txt",head=T)
sample_names<-read.table("sample_names_1559.txt",head=F)
names(pca_input)<-as.character(sample_names$V1)

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
d<-dist(new_pca2)
plot(hclust(d,"ward.D2"))
cut_clust<-cutree(hclust(d,"ward.D2"),4)


d<-dist(new_pca2)
pdf("Hierarhical_clustering_tree.pdf",width = 14,height = 7)
par(mfrow=c(1,1),mar=c(2,5,3,2)+0.1)
plot(hclust(d,"ward.D2"),xlab="",main="",cex=0.2,lwd=0.5)
abline(h=600, col="gray30",lty=3)
abline(h=330, col="gray30",lty=2)
text(c(710,1020,1420),c(620,620,620),c("Cluster 3","Cluster 2", "Cluster 1"),cex=0.7,col=c(2,3,4))
text(c(400,870,1020,1420),c(350,350,350),c("Cluster 3","Cluster 4","Cluster 2", "Cluster 1"),cex=0.7,col=c(2,5,3,4))
legend("topleft",col="grey30",lty=c(3,2),c("Cut with 3 branches", "Cut with 4 branches"),bty="n",cex=0.7)
dev.off()
cut_clust_4<-1+cutree(hclust(d,"ward.D2"),4)
cut_clust_3<-1+cutree(hclust(d,"ward.D2"),3)
pdf("Cluster_structure.pdf")
par(mfrow=c(2,2))
plot(sample_pca$x[,1],jitter(sample_pca$x[,2]),pch=19,cex=0.5, xlab="PC 1 (10%)",ylab="PC 2 (7%)",las=1,col=cut_clust_3,main="")
legend("topleft",pch=19,bty="n", c("Cluster 1", "Cluster 2", "Cluster 3"),col=c(4,3,2),cex=0.5)
mtext(side=2,at=40, line=2.7,"A",las=1,cex=2)
plot(sample_pca$x[,1],jitter(sample_pca$x[,2]),pch=19,cex=0.5, xlab="PC 1 (10%)",ylab="PC 2 (7%)",las=1,col=cut_clust_4,main="")
legend("topleft",pch=19,bty="n", c("Cluster 1", "Cluster 2", "Cluster 3","Cluster 4"),col=c(4,3,2,5),cex=0.5)
mtext(side=2,at=40, line=2.7,"B",las=1,cex=2)
plot(sample_pca$x[,1],jitter(sample_pca$x[,2]),pch=19,cex=0.5, xlab="PC 1 (10%)",ylab="PC 2 (7%)",las=1,main="",col=gray(missingness))
legend("topleft",col=c(gray(0.9),gray(0.1)),pch=19,bty="n",c("High missingness","Low missingness"),cex=0.5)
mtext(side=2,at=40, line=2.7,"C",las=1,cex=2)
plot(sample_pca$x[,1],jitter(sample_pca$x[,2]),pch=19,cex=0.5, xlab="PC 1 (10%)",ylab="PC 2 (7%)",las=1,main="",col=5)
points(sample_pca$x[sample_pca$x[,2] > 5,1],jitter(sample_pca$x[sample_pca$x[,2] > 5,2]),col=4,pch=19,cex=0.5)
points(sample_pca$x[sample_pca$x[,2] < -5 & sample_pca$x[,1] > 0,1],jitter(sample_pca$x[sample_pca$x[,2] < -5 & sample_pca$x[,1] > 0,2]),col=3,pch=19,cex=0.5)
points(sample_pca$x[sample_pca$x[,1] < -2,1],jitter(sample_pca$x[sample_pca$x[,1] < -2,2]),col=2,pch=19,cex=0.5)
legend("topleft",pch=19,bty="n", c("Cluster 1", "Cluster 2", "Cluster 3","Removed from further analysis"),col=c(4,3,2,5),cex=0.5)
#lines(y=c(5,5),x=c(-20,30),lty=3)
#lines(y=c(5,30),x=c(2,2),lty=3)
#lines(y=c(-5,-5),x=c(0,30),lty=3)
#lines(y=c(-5,-20),x=c(0,0),lty=3)
#lines(y=c(5,5),x=c(-20,-2),lty=3)
#lines(y=c(5,-25),x=c(-2,-2),lty=3)
mtext(side=2,at=40, line=2.7,"D",las=1,cex=2)
dev.off()


missingness<-apply(t(pca_input),1,function(x) length(which(x=="-")))
missingness<-missingness/14794

c(1:1559)[sample_pca$x[,2]>5] -> cluster_1
c(1:1559)[sample_pca$x[,2] < -5 & sample_pca$x[,1] > 0] -> cluster_2
c(1:1559)[sample_pca$x[,1] < -2] -> cluster_3
       
save.image("create_PCA.RData")

PCA_input_C1<-pca_input[,cluster_1]
PCA_input_C2<-pca_input[,cluster_2]
PCA_input_C3<-pca_input[,cluster_3]

PCA_C1<-rep(NA,dim(PCA_input_C1)[2]*dim(PCA_input_C1)[1])
dim(PCA_C1)<-dim(PCA_input_C1)
for (i in 1:dim(PCA_input_C1)[2]){
  PCA_C1[PCA_input_C1[,i]=="b",i]<-1
  PCA_C1[PCA_input_C1[,i]=="a",i]<-0
}
PCA_C1<-as.data.frame(t(PCA_C1))
apply(PCA_C1,2,mean,na.rm=T)->p.mean_C1
for(i in 1:dim(PCA_C1)[2]){
  PCA_C1[is.na(PCA_C1[,i]),i]<-p.mean_C1[i]
}
prcomp(~.,data=PCA_C1,na.action="na.exclude")->sample_pca_C1

PCA_C2<-rep(NA,dim(PCA_input_C2)[2]*dim(PCA_input_C2)[1])
dim(PCA_C2)<-dim(PCA_input_C2)
for (i in 1:dim(PCA_input_C2)[2]){
  PCA_C2[PCA_input_C2[,i]=="b",i]<-1
  PCA_C2[PCA_input_C2[,i]=="a",i]<-0
}
PCA_C2<-as.data.frame(t(PCA_C2))
apply(PCA_C2,2,mean,na.rm=T)->p.mean_C2
for(i in 1:dim(PCA_C2)[2]){
  PCA_C2[is.na(PCA_C2[,i]),i]<-p.mean_C2[i]
}
prcomp(~.,data=PCA_C2,na.action="na.exclude")->sample_pca_C2

PCA_C3<-rep(NA,dim(PCA_input_C3)[2]*dim(PCA_input_C3)[1])
dim(PCA_C3)<-dim(PCA_input_C3)
for (i in 1:dim(PCA_input_C3)[2]){
  PCA_C3[PCA_input_C3[,i]=="b",i]<-1
  PCA_C3[PCA_input_C3[,i]=="a",i]<-0
}
PCA_C3<-as.data.frame(t(PCA_C3))
apply(PCA_C3,2,mean,na.rm=T)->p.mean_C3
for(i in 1:dim(PCA_C3)[2]){
  PCA_C3[is.na(PCA_C3[,i]),i]<-p.mean_C3[i]
}
prcomp(~.,data=PCA_C3,na.action="na.exclude")->sample_pca_C3

summary(sample_pca_C1)$importance[,1:5]
summary(sample_pca_C2)$importance[,1:5]
summary(sample_pca_C3)$importance[,1:5]

pdf("Per_cluster_PCA.pdf",height = 7,width = 2)
  par(mfrow=c(3,1))
  plot(sample_pca_C1$x[,1],jitter(sample_pca_C1$x[,2]),pch=19,cex=0.5, xlab="PC 1 (3.6%)",ylab="PC 2 (3.3%)",las=1,main="Cluster 1",col=1)
  plot(sample_pca_C2$x[,1],jitter(sample_pca_C2$x[,2]),pch=19,cex=0.5, xlab="PC 1 (3.8%)",ylab="PC 2 (3.7%)",las=1,main="Cluster 2",col=1)
  plot(sample_pca_C3$x[,1],jitter(sample_pca_C3$x[,2]),pch=19,cex=0.5, xlab="PC 1 (3.9%)",ylab="PC 2 (3.5%)",las=1,main="Cluster 3",col=1)
dev.off()

save.image("create_PCA.RData")