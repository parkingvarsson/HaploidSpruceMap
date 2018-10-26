##Import data
read.table("C1_all_marker_maps.txt",head=T)->C1_full
read.table("C2_all_marker_maps.txt",head=T)->C2_full
read.table("C3_all_marker_maps.txt",head=T)->C3_full

C1_full$scaffold<- as.factor(t(as.data.frame(strsplit(as.character(C1_full$Marker),split=":")))[,1])
C1_full$probe<- t(as.data.frame(strsplit(as.character(C1_full$Marker),split=":")))[,2]
C2_full$scaffold<- as.factor(t(as.data.frame(strsplit(as.character(C2_full$Marker),split=":")))[,1])
C2_full$probe<- t(as.data.frame(strsplit(as.character(C2_full$Marker),split=":")))[,2]
C3_full$scaffold<- as.factor(t(as.data.frame(strsplit(as.character(C3_full$Marker),split=":")))[,1])
C3_full$probe<- t(as.data.frame(strsplit(as.character(C3_full$Marker),split=":")))[,2]

C1_scaf_vs_LG<-table(C1_full$scaffold,C1_full$LG)
C2_scaf_vs_LG<-table(C2_full$scaffold,C2_full$LG)
C3_scaf_vs_LG<-table(C3_full$scaffold,C3_full$LG)

read.table("C1_probe-pos_info.txt",head=T)->C1_probePos
read.table("C2_probe-pos_info.txt",head=T)->C2_probePos
read.table("C3_probe-pos_info.txt",head=T)->C3_probePos

C1_full<-merge(C1_full, C1_probePos[,2:3],by.x="Marker",by.y="scaf_probe")
C1_full<-C1_full[order(C1_full$LG,C1_full$Position),]
C2_full<-merge(C2_full, C2_probePos[,2:3],by.x="Marker",by.y="scaf_probe")
C2_full<-C2_full[order(C2_full$LG,C2_full$Position),]
C3_full<-merge(C3_full, C3_probePos[,2:3],by.x="Marker",by.y="scaf_probe")
C3_full<-C3_full[order(C3_full$LG,C3_full$Position),]


C1_full$rev_position<-NA
C2_full$rev_position<-NA
C3_full$rev_position<-NA
for(i in 1:12){
  C1_full$rev_position[C1_full$LG==i]<-max(C1_full$Position[C1_full$LG==i])-C1_full$Position[C1_full$LG==i]
  C2_full$rev_position[C2_full$LG==i]<-max(C2_full$Position[C2_full$LG==i])-C2_full$Position[C2_full$LG==i]
  C3_full$rev_position[C3_full$LG==i]<-max(C3_full$Position[C3_full$LG==i])-C3_full$Position[C3_full$LG==i]
}


InterSplit_scaf_C1<-rownames(C1_scaf_vs_LG[which(apply(C1_scaf_vs_LG,1,function(x) length(which(x > 0))) >1),])
InterSplit_scaf_C2<-rownames(C2_scaf_vs_LG[which(apply(C2_scaf_vs_LG,1,function(x) length(which(x > 0))) >1),])
InterSplit_scaf_C3<-rownames(C3_scaf_vs_LG[which(apply(C3_scaf_vs_LG,1,function(x) length(which(x > 0))) >1),])


#Which LGs are correlated between clusters
comp_C1<-NULL
comp_C2<-NULL
comp_C3<-NULL
comp_C1_C2<-NULL
comp_C1_C3<-NULL
comp_C2_C3<-NULL
k<-1
for(i in 1:12){
  for(j in 1:12){
    comp_C1[k]<-length(as.character(unique(C1_full$scaffold[C1_full$LG==i & C1_full$scaffold %in% C1_full$scaffold[C1_full$LG==j]])))
    comp_C2[k]<-length(as.character(unique(C2_full$scaffold[C2_full$LG==i & C2_full$scaffold %in% C2_full$scaffold[C2_full$LG==j]])))
    comp_C3[k]<-length(as.character(unique(C3_full$scaffold[C3_full$LG==i & C3_full$scaffold %in% C3_full$scaffold[C3_full$LG==j]])))
    comp_C1_C2[k]<-length(as.character(unique(C1_full$scaffold[C1_full$LG==i & C1_full$scaffold %in% C2_full$scaffold[C2_full$LG==j]])))
    comp_C1_C3[k]<-length(as.character(unique(C1_full$scaffold[C1_full$LG==i & C1_full$scaffold %in% C3_full$scaffold[C3_full$LG==j]])))
    comp_C2_C3[k]<-length(as.character(unique(C2_full$scaffold[C2_full$LG==i & C2_full$scaffold %in% C3_full$scaffold[C3_full$LG==j]])))
    
    k<-k+1
  }  
}
comp<-matrix(NA,36,36)
comp[1:12,1:12]<-matrix(comp_C1,12,12)
comp[13:24,13:24]<-matrix(comp_C2,12,12)
comp[25:36,25:36]<-matrix(comp_C3,12,12)
comp[13:24,1:12]<-matrix(comp_C1_C2,12,12)
comp[25:36,1:12]<-matrix(comp_C1_C3,12,12)
comp[25:36,13:24]<-matrix(comp_C2_C3,12,12)
comp[upper.tri(comp,diag = F)]<-NA

## Plot heatmap of cluster comparison
install.packages("pheatmap")
library(pheatmap)
m<-comp
labels<-c("C1_1","C1_2","C1_3","C1_4","C1_5","C1_6","C1_7","C1_8","C1_9","C1_0","C1_11","C1_12","C2_1","C2_2","C2_3","C2_4","C2_5","C2_6","C2_7","C2_8","C2_9","C2_10","C2_11","C2_12","C3_1","C3_2","C3_3","C3_4","C3_5","C3_6","C3_7","C3_8","C3_9","C3_10","C3_11","C3_12")
pdf("Scaffolds_shared_between_cluster.pdf")
pheatmap(comp,cluster_rows=FALSE,cluster_cols=FALSE,main = "Shared scaffolds between Cluster LGs",gaps_row = c(12,24),gaps_col = c(12,24),display_numbers = matrix(ifelse(is.na(m)|m<1,"",m),nrow(m)),na_col = "white",labels_col = labels, labels_row = labels,fontsize_number = 6)
dev.off()

### Table with LG correlation
LG_comp<-NULL
LG_comp$LG<-c(1,2,3,4,5,6,7,8,9,10,11,12)
LG_comp$C1<-c(1,2,3,4,5,6,7,8,9,10,11,12)
LG_comp$C2<-c(11,12,1,3,9,5,6,8,4,10,2,7)
LG_comp$C3<-c(10,1,2,5,3,7,8,11,6,12,4,9)
LG_comp<-as.data.frame(LG_comp)



###Merge maps into the same data frame
ClusterMaps_merged<-merge(C1_full,C2_full,by=c("Marker","scaffold","probe"),all = T)
ClusterMaps_merged<-merge(ClusterMaps_merged,C3_full,by=c("Marker","scaffold","probe"),all = T)
names(ClusterMaps_merged)<-c("Marker",  "scaffold", "probe", "LG.C1",  "gen_pos.C1",  "phys_pos.C1","revGen_pos_C1", "LG.C2",  "gen_pos.C2", "phys_pos.C2","revGen_pos_C2", "LG.C3", "gen_pos.C3", "phys_pos.C3","revGen_pos_C3")
ClusterMaps_merged$LG<-NA
for(i in 1:dim(ClusterMaps_merged)[1]){
  if(!is.na(ClusterMaps_merged$LG.C1[i])){
    ClusterMaps_merged$LG[i]<-LG_comp$LG[LG_comp$C1==ClusterMaps_merged$LG.C1[i]]
  }else{
    if(!is.na(ClusterMaps_merged$LG.C2[i])){
      ClusterMaps_merged$LG[i]<-LG_comp$LG[LG_comp$C2==ClusterMaps_merged$LG.C2[i]]
    } else{
      ClusterMaps_merged$LG[i]<-LG_comp$LG[LG_comp$C3==ClusterMaps_merged$LG.C3[i]]
    }
  }  
}
ClusterMaps_merged<-ClusterMaps_merged[order(ClusterMaps_merged$LG,ClusterMaps_merged$gen_pos.C1),c(1:3,16,4:15)]

ClusterMaps_merged$GenPos_C1<-ClusterMaps_merged$gen_pos.C1

ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==1]<-ClusterMaps_merged$revGen_pos_C2[ClusterMaps_merged$LG==1]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==2]<-ClusterMaps_merged$revGen_pos_C2[ClusterMaps_merged$LG==2]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==3]<-ClusterMaps_merged$gen_pos.C2[ClusterMaps_merged$LG==3]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==4]<-ClusterMaps_merged$revGen_pos_C2[ClusterMaps_merged$LG==4]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==5]<-ClusterMaps_merged$gen_pos.C2[ClusterMaps_merged$LG==5]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==6]<-ClusterMaps_merged$gen_pos.C2[ClusterMaps_merged$LG==6]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==7]<-ClusterMaps_merged$revGen_pos_C2[ClusterMaps_merged$LG==7]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==8]<-ClusterMaps_merged$gen_pos.C2[ClusterMaps_merged$LG==8]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==9]<-ClusterMaps_merged$revGen_pos_C2[ClusterMaps_merged$LG==9]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==10]<-ClusterMaps_merged$revGen_pos_C2[ClusterMaps_merged$LG==10]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==11]<-ClusterMaps_merged$revGen_pos_C2[ClusterMaps_merged$LG==11]
ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==12]<-ClusterMaps_merged$gen_pos.C2[ClusterMaps_merged$LG==12]

ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==1]<-ClusterMaps_merged$revGen_pos_C3[ClusterMaps_merged$LG==1]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==2]<-ClusterMaps_merged$revGen_pos_C3[ClusterMaps_merged$LG==2]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==3]<-ClusterMaps_merged$gen_pos.C3[ClusterMaps_merged$LG==3]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==4]<-ClusterMaps_merged$gen_pos.C3[ClusterMaps_merged$LG==4]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==5]<-ClusterMaps_merged$revGen_pos_C3[ClusterMaps_merged$LG==5]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==6]<-ClusterMaps_merged$gen_pos.C3[ClusterMaps_merged$LG==6]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==7]<-ClusterMaps_merged$revGen_pos_C3[ClusterMaps_merged$LG==7]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==8]<-ClusterMaps_merged$gen_pos.C3[ClusterMaps_merged$LG==8]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==9]<-ClusterMaps_merged$revGen_pos_C3[ClusterMaps_merged$LG==9]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==10]<-ClusterMaps_merged$gen_pos.C3[ClusterMaps_merged$LG==10]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==11]<-ClusterMaps_merged$revGen_pos_C3[ClusterMaps_merged$LG==11]
ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==12]<-ClusterMaps_merged$revGen_pos_C3[ClusterMaps_merged$LG==12]

ClusterMaps_merged<-ClusterMaps_merged[,c(1:5,17,7,9,18,11,13,19,15)]

par(mfrow=c(3,4))
for(i in 1:12){
  plot(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==i],ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==i],pch=19,col="black",cex=0.5,xlab="First cluster Position (cM)",ylab="Second cluster Position (cM)",main=paste("LG",i),las=2,xlim=c(0,440),ylim=c(0,440))
  points(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==i],ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==i],pch=19,col="red",cex=0.5)
  points(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==i],ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==i],pch=19,col="blue",cex=0.5)
  legend("topleft",c("C1 - C2", "C1 - C3", "C3 - C2"),col=c("black","red","blue"),cex=0.7,pch=19,bty="n")  
}


###Interspit scaffolds

Intersplit_scaf_maps<-ClusterMaps_merged[ClusterMaps_merged$scaffold %in% unique(c(InterSplit_scaf_C1,InterSplit_scaf_C2,InterSplit_scaf_C3)),]
Intersplit_scaf_maps<-Intersplit_scaf_maps[order(Intersplit_scaf_maps$scaffold),]


###Intrasplit scaffolds

IntraLG_scafPos<-NULL
for(i in 1:length(levels(ClusterMaps_merged$scaffold))){
  df<-ClusterMaps_merged[ClusterMaps_merged$scaffold==levels(ClusterMaps_merged$scaffold)[i],]
  for(j in 1:length(unique(df$LG))){
    if(dim(df[df$LG==unique(df$LG)[j],])[1]>1){
      tmp<-df[df$LG==unique(df$LG)[j],]
      tmp<-tmp[order(tmp$probe),]
      for(k in 1:(dim(tmp)[1]-1)){
        for(l in 2:dim(tmp)[1]){
          if(l >k){
            IntraLG_scafPos<-rbind(IntraLG_scafPos,cbind(unique(as.character(tmp$scaffold)),unique(tmp$LG),paste(tmp$probe[k],"_",tmp$probe[l],sep=""),abs(tmp$GenPos_C1[k]-tmp$GenPos_C1[l]),
                                                   abs(tmp$phys_pos.C1[k]-tmp$phys_pos.C1[l]),abs(tmp$GenPos_C2[k]-tmp$GenPos_C2[l]),
                                                   abs(tmp$phys_pos.C2[k]-tmp$phys_pos.C2[l]),abs(tmp$GenPos_C3[k]-tmp$GenPos_C3[l]),
                                                   abs(tmp$phys_pos.C3[k]-tmp$phys_pos.C3[l])))
          }
        }
      }  
    }  
  }
}
rm(df,tmp)

IntraLG_scafPos<-as.data.frame(IntraLG_scafPos)
names(IntraLG_scafPos)<-c("Scaffold","LG","Probe_comp","C1_GenDist","C1_PhysDist","C2_GenDist","C2_PhysDist","C3_GenDist","C3_PhysDist")
IntraLG_scafPos$LG<-as.integer(as.character(IntraLG_scafPos$LG))
IntraLG_scafPos$C1_GenDist<-as.numeric(as.character(IntraLG_scafPos$C1_GenDist))
IntraLG_scafPos$C2_GenDist<-as.numeric(as.character(IntraLG_scafPos$C2_GenDist))
IntraLG_scafPos$C3_GenDist<-as.numeric(as.character(IntraLG_scafPos$C3_GenDist))
IntraLG_scafPos$C1_PhysDist<-as.numeric(as.character(IntraLG_scafPos$C1_PhysDist))
IntraLG_scafPos$C2_PhysDist<-as.numeric(as.character(IntraLG_scafPos$C2_PhysDist))
IntraLG_scafPos$C3_PhysDist<-as.numeric(as.character(IntraLG_scafPos$C3_PhysDist))

IntraLG_scafPos$C1_GenToPhys<-IntraLG_scafPos$C1_GenDist/IntraLG_scafPos$C1_PhysDist
IntraLG_scafPos$C2_GenToPhys<-IntraLG_scafPos$C2_GenDist/IntraLG_scafPos$C2_PhysDist
IntraLG_scafPos$C3_GenToPhys<-IntraLG_scafPos$C3_GenDist/IntraLG_scafPos$C3_PhysDist
IntraLG_scafPos<-IntraLG_scafPos[which(apply(IntraLG_scafPos,1,function(x) sum(is.na(x)))!=9),]

save.image("Cluster_map_comparisons.RData")

### plot 
par(mfrow=c(3,2),mar=c(5, 4, 4, 2) + 0.1)
plot(IntraLG_scafPos$C1_PhysDist,log10(IntraLG_scafPos$C1_GenDist),main="C1 - Genetic distance vs Physical distance",xlab="Physical distance (bp)",ylab= expression(paste(log[10],"(Genetic distance (cM))")),las=1,ylim=c(-5,3),pch=19,col=rgb(0.4,0.4,0.4,alpha=0.5))
plot(density(IntraLG_scafPos$C1_GenToPhys,bw=0.001,na.rm=T),main="C1 - Genetic distance/Physical distance",xlim=c(-0.01,0.05))
abline(v=quantile(IntraLG_scafPos$C1_GenToPhys,0.95,na.rm=T)[[1]],lty=3,col="red")
plot(IntraLG_scafPos$C2_PhysDist,log10(IntraLG_scafPos$C2_GenDist),main="C2 - Genetic distance vs Physical distance",xlab="Physical distance (bp)",ylab= expression(paste(log[10],"(Genetic distance (cM))")),las=1,ylim=c(-5,3),pch=19,col=rgb(0.4,0.4,0.4,alpha=0.5))
plot(density(IntraLG_scafPos$C2_GenToPhys,bw=0.001,na.rm=T),main="C2 - Genetic distance/Physical distance",xlim=c(-0.01,0.05))
abline(v=quantile(IntraLG_scafPos$C2_GenToPhys,0.95,na.rm=T)[[1]],lty=3,col="red")
plot(IntraLG_scafPos$C3_PhysDist,log10(IntraLG_scafPos$C3_GenDist),main="C3 - Genetic distance vs Physical distance",xlab="Physical distance (bp)",ylab=expression(paste(log[10],"(Genetic distance (cM))")),las=1,ylim=c(-5,3),pch=19,col=rgb(0.4,0.4,0.4,alpha=0.5))
plot(density(IntraLG_scafPos$C3_GenToPhys,bw=0.001,na.rm=T),main="C3 - Genetic distance/Physical distance",xlim=c(-0.01,0.05))
abline(v=quantile(IntraLG_scafPos$C3_GenToPhys,0.95,na.rm=T)[[1]],lty=3,col="red")


par(mfcol=c(3,3),mar=c(5, 4, 4, 2) + 0.1)
plot(IntraLG_scafPos$C1_GenDist,IntraLG_scafPos$C2_GenDist,col=rgb(0.4,0.4,0.4,alpha=0.6),pch=19,xlab="C1",ylab = "C2",main="Genetic distance comparison between C1 and C2",las=1)#,xlim=c(0,20),ylim=c(0,20))
plot(IntraLG_scafPos$C1_GenDist,IntraLG_scafPos$C3_GenDist,col=rgb(0.4,0.4,0.4,alpha=0.6),pch=19,xlab="C1",ylab = "C3",main="Genetic distance comparison between C1 and C3",las=1)#,xlim=c(0,20),ylim=c(0,20))
plot(IntraLG_scafPos$C2_GenDist,IntraLG_scafPos$C3_GenDist,col=rgb(0.4,0.4,0.4,alpha=0.6),pch=19,xlab="C2",ylab = "C3",main="Genetic distance comparison between C2 and C3",las=1)#,xlim=c(0,20),ylim=c(0,20))
plot(IntraLG_scafPos$C1_GenDist,IntraLG_scafPos$C2_GenDist,col=rgb(0.4,0.4,0.4,alpha=0.6),pch=19,xlab="C1",ylab = "C2",main="Genetic distance comparison between C1 and C2",las=1,xlim=c(0,20),ylim=c(0,20))
lines(c(3,20),c(3,3),col="red",lty=3)
lines(c(3,3),c(3,20),col="red",lty=3)
plot(IntraLG_scafPos$C1_GenDist,IntraLG_scafPos$C3_GenDist,col=rgb(0.4,0.4,0.4,alpha=0.6),pch=19,xlab="C1",ylab = "C3",main="Genetic distance comparison between C1 and C3",las=1,xlim=c(0,20),ylim=c(0,20))
lines(c(3,20),c(3,3),col="red",lty=3)
lines(c(3,3),c(3,20),col="red",lty=3)
plot(IntraLG_scafPos$C2_GenDist,IntraLG_scafPos$C3_GenDist,col=rgb(0.4,0.4,0.4,alpha=0.6),pch=19,xlab="C2",ylab = "C3",main="Genetic distance comparison between C2 and C3",las=1,xlim=c(0,20),ylim=c(0,20))
lines(c(3,20),c(3,3),col="red",lty=3)
lines(c(3,3),c(3,20),col="red",lty=3)
scatterplot3d(IntraLG_scafPos$C1_GenDist,IntraLG_scafPos$C2_GenDist,IntraLG_scafPos$C3_GenDist,color = rgb(0.4,0.4,0.4,alpha=0.6),pch=19,xlab="C1",ylab="C2",zlab="C3",las=1)#,xlim=c(0,20),ylim=c(0,20),zlim=c(0,20))
scatterplot3d(IntraLG_scafPos$C1_GenDist,IntraLG_scafPos$C2_GenDist,IntraLG_scafPos$C3_GenDist,color = rgb(0.4,0.4,0.4,alpha=0.6),pch=19,xlab="C1",ylab="C2",zlab="C3",las=1,xlim=c(0,20),ylim=c(0,20),zlim=c(0,20))



gap_density<-NULL
gap_density$C1<-hist(IntraLG_scafPos$C1_GenDist,breaks=140)$density
gap_density$C2<-hist(IntraLG_scafPos$C2_GenDist,breaks=220)$density
gap_density$C3<-hist(IntraLG_scafPos$C3_GenDist,breaks=180)$density

plot.new()
plot.window(xlim=c(1,20),ylim=c(0.8,1))
abline(v=seq(1,20,1),lty=3,col="grey80")
abline(h=seq(0.8,1,0.01),lty=6,col="grey80")
points(df$C2,df$quantile,pch=19,cex=0.5,col=rgb(0,1,0,alpha=0.8),type="l")
points(df$C1,df$quantile,pch=19,cex=0.5,col=rgb(1,0,0,alpha=0.8),type="l")
points(df$C3,df$quantile,pch=19,cex=0.5,col=rgb(0,0,1,alpha=0.8),type="l")
axis(1);axis(2,las=1);box()
mtext(side=1,text="Probe-pair gap (cM)",line=2.5)
mtext(side=2,text="Fraction of probe-pairs",line=3)

####For LPmerge
LGI<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==1 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==1 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==1 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==1 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==1 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==1 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGII<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==2 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==2 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==2 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==2 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==2 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==2 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGIII<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==3 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==3 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==3 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==3 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==3 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==3 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGIV<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==4 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==4 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==4 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==4 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==4 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==4 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGV<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==5 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==5 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==5 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==5 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==5 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==5 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGVI<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==6 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==6 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==6 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==6 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==6 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==6 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGVII<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==7 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==7 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==7 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==7 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==7 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==7 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGVIII<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==8 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==8 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==8 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==8 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==8 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==8 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGIX<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==9 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==9 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==9 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==9 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==9 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==9 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGX<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==10 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==10 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==10 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==10 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==10 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==10 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGXI<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==11 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==11 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==11 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==11 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==11 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==11 & !is.na(ClusterMaps_merged$GenPos_C1)]),])
LGXII<-list(C3=ClusterMaps_merged[ClusterMaps_merged$LG==12 & !is.na(ClusterMaps_merged$GenPos_C3),c(1,12)][order(ClusterMaps_merged$GenPos_C3[ClusterMaps_merged$LG==12 & !is.na(ClusterMaps_merged$GenPos_C3)]),],C2=ClusterMaps_merged[ClusterMaps_merged$LG==12 & !is.na(ClusterMaps_merged$GenPos_C2),c(1,9)][order(ClusterMaps_merged$GenPos_C2[ClusterMaps_merged$LG==12 & !is.na(ClusterMaps_merged$GenPos_C2)]),],C1=ClusterMaps_merged[ClusterMaps_merged$LG==12 & !is.na(ClusterMaps_merged$GenPos_C1),c(1,6)][order(ClusterMaps_merged$GenPos_C1[ClusterMaps_merged$LG==12 & !is.na(ClusterMaps_merged$GenPos_C1)]),])

save(LGI,LGII,LGIII,LGIV,LGV,LGVI,LGVII,LGVIII,LGIX,LGX,LGXI,LGXII,file="LGs_for_consensus.RData")

C3_gaps<-NULL
C3_gaps$LG<-NULL
C3_gaps$gap<-NULL
for(i in 1:12){
for(j in 2:dim(C3_full[C3_full$LG==i,])[1]){
  df<-C3_full[C3_full$LG==i,]
  C3_gaps$LG<-c(C3_gaps$LG,i)
  C3_gaps$gap<-c(C3_gaps$gap,df$Position[j]-df$Position[j-1])
  }
}  

C3_gaps<-as.data.frame(C3_gaps)

for(i in 1:12){
  print(length(which(C3_gaps$gap[C3_gaps$LG==i] < 5.01))/length(C3_gaps$gap[C3_gaps$LG==i])*100)
}

C2_gaps<-NULL
C2_gaps$LG<-NULL
C2_gaps$gap<-NULL
for(i in 1:12){
  for(j in 2:dim(C2_full[C2_full$LG==i,])[1]){
    df<-C2_full[C2_full$LG==i,]
    C2_gaps$LG<-c(C2_gaps$LG,i)
    C2_gaps$gap<-c(C2_gaps$gap,df$Position[j]-df$Position[j-1])
  }
}  

C2_gaps<-as.data.frame(C2_gaps)

for(i in 1:12){
  print(length(which(C2_gaps$gap[C2_gaps$LG==i] < 5.01))/length(C2_gaps$gap[C2_gaps$LG==i])*100)
}

C1_gaps<-NULL
C1_gaps$LG<-NULL
C1_gaps$gap<-NULL
for(i in 1:12){
  for(j in 2:dim(C1_full[C1_full$LG==i,])[1]){
    df<-C1_full[C1_full$LG==i,]
    C1_gaps$LG<-c(C1_gaps$LG,i)
    C1_gaps$gap<-c(C1_gaps$gap,df$Position[j]-df$Position[j-1])
  }
}  

C1_gaps<-as.data.frame(C1_gaps)

for(i in 1:12){
  print(length(C1_gaps$gap[C1_gaps$LG==i & C1_gaps$gap <5.01])/length(C1_gaps$gap[C1_gaps$LG==i])*100)
}
