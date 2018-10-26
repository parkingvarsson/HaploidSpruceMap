load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_1/C1_cleaned.RData")

library(BatchMap)
library(onemap)

load("C1_LG_1_rec16_10wRippled_Map.RData")
C1_LG1<-map$Map
pdf("C1_LG1_heatmap.pdf")
  rf_graph_table(C1_LG1,scale=2,inter=F,main="C1_LG1",axis.cex = 0.4)
dev.off()
load("C1_LG_2_rec16_10wRippled_Map.RData")
C1_LG2<-map$Map
pdf("C1_LG2_heatmap.pdf")
rf_graph_table(C1_LG2,scale=2,inter=F,main="C1_LG2",axis.cex = 0.4)
dev.off()
load("C1_LG_3_rec16_10wRippled_Map.RData")
C1_LG3<-map$Map
pdf("C1_LG3_heatmap.pdf")
rf_graph_table(C1_LG3,scale=2,inter=F,main="C1_LG3",axis.cex = 0.4)
dev.off()
load("C1_LG_4_rec16_10wRippled_Map.RData")
C1_LG4<-map$Map
pdf("C1_LG4_heatmap.pdf")
rf_graph_table(C1_LG4,scale=2,inter=F,main="C1_LG4",axis.cex = 0.4)
dev.off()
load("C1_LG_5_rec16_10wRippled_Map.RData")
C1_LG5<-map$Map
pdf("C1_LG5_heatmap.pdf")
rf_graph_table(C1_LG5,scale=2,inter=F,main="C1_LG5",axis.cex = 0.4)
dev.off()
load("C1_LG_6_rec16_10wRippled_Map.RData")
C1_LG6<-map$Map
pdf("C1_LG6_heatmap.pdf")
rf_graph_table(C1_LG6,scale=2,inter=F,main="C1_LG6",axis.cex = 0.4)
dev.off()
load("C1_LG_7_rec16_10wRippled_Map.RData")
C1_LG7<-map$Map
pdf("C1_LG7_heatmap.pdf")
rf_graph_table(C1_LG7,scale=2,inter=F,main="C1_LG7",axis.cex = 0.4)
dev.off()
load("C1_LG_8_rec16_10wRippled_Map.RData")
C1_LG8<-map$Map
pdf("C1_LG8_heatmap.pdf")
rf_graph_table(C1_LG8,scale=2,inter=F,main="C1_LG8",axis.cex = 0.4)
dev.off()
load("C1_LG_9_rec16_10wRippled_Map.RData")
C1_LG9<-map$Map
pdf("C1_LG9_heatmap.pdf")
rf_graph_table(C1_LG9,scale=2,inter=F,main="C1_LG9",axis.cex = 0.4)
dev.off()
load("C1_LG_10_rec16_10wRippled_Map.RData")
C1_LG10<-map$Map
pdf("C1_LG10_heatmap.pdf")
rf_graph_table(C1_LG10,scale=2,inter=F,main="C1_LG10",axis.cex = 0.4)
dev.off()
load("C1_LG_11_rec16_10wRippled_Map.RData")
C1_LG11<-map$Map
pdf("C1_LG11_heatmap.pdf")
rf_graph_table(C1_LG11,scale=2,inter=F,main="C1_LG11",axis.cex = 0.4)
dev.off()
load("C1_LG_12_rec16_10wRippled_Map.RData")
C1_LG12<-map$Map
pdf("C1_LG12_heatmap.pdf")
rf_graph_table(C1_LG12,scale=2,inter=F,main="C1_LG12",axis.cex = 0.4)
dev.off()

save.image("C1_cleaned_maps.RData")

### WRITE MAPS TO FILE
write.map(list(C1_LG1,C1_LG2,C1_LG3,C1_LG4,C1_LG5,C1_LG6,C1_LG7,C1_LG8,C1_LG9,C1_LG10,C1_LG11,C1_LG12),"C1_bin_maps.txt")

bin_maps<-read.table("C1_bin_maps.txt",head=F)
names(bin_maps)<-c("LG","Marker","Position")


segregation_clean<-test_segregation(C1.out_clean)
plot.onemap_segreg_test(segregation_clean)
pdf("Segregation_distortion_C1.pdf")
  plot_segreg_test(segregation_clean,Title="Cluster 1")
dev.off()
segregation_clean_table<-print.onemap.segreg.test(segregation_clean)
bin_map_with_segregation_pattern<-merge(bin_maps,segregation_clean_table,by.x="Marker",by.y="Marker")
bin_map_with_segregation_pattern<-bin_map_with_segregation_pattern[order(bin_map_with_segregation_pattern$LG,bin_map_with_segregation_pattern$Position),]

par(mfrow=c(3,4))
for(i in 1:12){
  plot(bin_map_with_segregation_pattern$Position[bin_map_with_segregation_pattern$LG==i],-log10(bin_map_with_segregation_pattern$`p-value`[bin_map_with_segregation_pattern$LG==i]),pch=19,cex=0.4,main=paste("C1 LG ",i),ylab=expression(paste(-log[10],"(p-value)")),xlab="Genetic position (cM)",las=2)
  abline(h=-log10(Bonferroni_alpha(segregation_clean)),col="red")
}



Full_maps<-NULL
for(i in 1:length(bins$bins)){
  Full_maps<-rbind(Full_maps,as.data.frame(cbind(rep(bin_maps$LG[bin_maps$Marker==names(bins$bins[i])],length(rownames(bins$bins[[i]]))),rownames(bins$bins[[i]]),rep(bin_maps$Position[bin_maps$Marker==names(bins$bins[i])],length(rownames(bins$bins[[i]]))))))
}
names(Full_maps)<-c("LG","Marker","Position")
Full_maps$LG<-as.integer(as.character(Full_maps$LG))
Full_maps$Position<-as.numeric(as.character(Full_maps$Position))

Full_maps<-Full_maps[order(Full_maps$LG,Full_maps$Position),]

write.table(Full_maps,"C1_all_marker_maps.txt",col.names = T,row.names = F,quote=F, sep="\t")

Full_maps$Scaffold<-as.factor(t(as.data.frame(strsplit(as.character(Full_maps$Marker),split = ":")))[,1])
Full_maps$Probe<-t(as.data.frame(strsplit(as.character(Full_maps$Marker),split = ":")))[,2]

Full_maps<-Full_maps[order("LG","Position"),]

C1_scaf_vs_C1_LG<-table(Full_maps$Scaffold,Full_maps$LG)

C1_scaf_vs_C1_LG[which(apply(C1_scaf_vs_C1_LG,1,function(x) length(which(x > 0)))>1),]
dim(C1_scaf_vs_C1_LG[which(apply(C1_scaf_vs_C1_LG,1,function(x) length(which(x > 0)))>1),])

IntraLG_scafPos<-NULL
for(i in 1:length(levels(Full_maps$Scaffold))){
  df<-Full_maps[Full_maps$Scaffold==levels(Full_maps$Scaffold)[i],]
  for(j in 1:length(unique(df$LG))){
    if(dim(df[df$LG==unique(df$LG)[j],])[1]>1){
        IntraLG_scafPos<-rbind(IntraLG_scafPos,cbind(unique(as.character(df$Scaffold)),unique(df$LG)[j],max(df$Position[df$LG==unique(df$LG)[j]])-min(df$Position[df$LG==unique(df$LG)[j]])))
    }
  }
}
IntraLG_scafPos<-as.data.frame(IntraLG_scafPos)
names(IntraLG_scafPos)<-c("Scaffold","LG","Distance")
IntraLG_scafPos$LG<-as.integer(as.character(IntraLG_scafPos$LG))
IntraLG_scafPos$Distance<-as.numeric(as.character(IntraLG_scafPos$Distance))

rm(df,i,j)
save.image("C1_cleaned_maps.RData")


### Marker distribution
par(mfrow=c(3,4))
for(i in 1:12){
  tmp<-density(Full_maps$Position[Full_maps$LG==i],from=min(Full_maps$Position[Full_maps$LG==i]),to=max(Full_maps$Position[Full_maps$LG==i]),bw=10)
  plot(tmp,main=paste("C1 LG ",i),ylab="Marker density", xlab="Genetic position (cM)")
  abline(h=mean(tmp$y), col="red")
  abline(h=mean(tmp$y)+sd(tmp$y),col="red", lty=3)
  abline(h=mean(tmp$y)-sd(tmp$y),col="red", lty=3)
}

###Largest gap in each LG
gap<-NULL
gap$LG<-seq(1,12,1)
gap$gap<-rep(0,12)
gap<-as.data.frame(gap)
for(i in 1:12){
  tmp<-Full_maps[Full_maps$LG==i,]
  for(j in 2:dim(tmp)[1]){
    if(tmp$Position[j]-tmp$Position[j-1] > gap$gap[gap$LG==i]){
      gap$gap[gap$LG==i]<-tmp$Position[j]-tmp$Position[j-1]
    }
  }
}