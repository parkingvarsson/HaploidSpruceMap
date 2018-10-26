#### Create Figure showing the number of markers per scaffold and how many of the scaffolds that are wrongly assembled
read.table("../Consensus_maps.txt",head=T)->maps
read.table("../../C1_probe-pos_info.txt",head=T)->C1_phys_pos
read.table("../../C2_probe-pos_info.txt",head=T)->C2_phys_pos
read.table("../../C3_probe-pos_info.txt",head=T)->C3_phys_pos

maps<-merge(maps,C3_phys_pos[,2:3],by.x="Marker",by.y = "scaf_probe",all = T)
maps<-merge(maps,C2_phys_pos[,2:3],by.x="Marker",by.y = "scaf_probe",all = T)
maps<-merge(maps,C1_phys_pos[,2:3],by.x="Marker",by.y = "scaf_probe",all = T)
names(maps)<-c("Marker","LG","Consensus","C3","C2","C1","Scaffold","Probe","C3_pos","C2_pos","C1_pos")

read.table("inter_split_scaffolds.txt",head=F)->inter_split
read.table("intra_split_scaffolds.txt",head=F)->intra_split

map_intra_split<- maps[maps$Scaffold %in% intra_split$V1,]
map_intra_split<-drop.levels(map_intra_split)
map_intra_split$Consensus_pos<- apply(map_intra_split[9:11],1,function(x) mean(x,na.rm=T))

read.table("Pabies1.0-all.phase.changed.gff3",head=F)->All_models
library(gdata)
All_models<-drop.levels(All_models[All_models$V1 %in% as.character(map_intra_split$Scaffold),])
read.table("haploid_scaffolds.txt",head=F)->scaf_length
read.table("intra_split_Scaffolds.gaps.bed",head=F)->gaps

pdf("IntraSplitScaffolds_Map_Comparisons.pdf")
par(mfrow=c(3,3))
plot.new()
legend("center",pch=c(19,17,18,6),lwd=c(2,1,1,1),col=c("black","green","red","blue"),lty=c(1,2,3,4),c("Consensus","C3","C2","C1"),bty="n",cex=1.5)
for(i in 1:length(levels(map_intra_split$Scaffold))){
  tmp<-map_intra_split[map_intra_split$Scaffold==levels(map_intra_split$Scaffold)[i],]
  tmp1<-All_models[as.character(All_models$V1)==as.character(unique(tmp$Scaffold)) & All_models$V3 %in% c("gene","exon"),c(1,3,4,5)]
  tmp2<-gaps[as.character(gaps$V1)==as.character(unique(tmp$Scaffold)),]
  plot(1,type="n", xlim=c(0,scaf_length$V2[scaf_length$V1==as.character(unique(tmp$Scaffold))]),ylim=c(min(tmp[,3:6],na.rm=T)-7,
          max(tmp[,3:6],na.rm=T)),pch=19,cex=0.5,col="black",
       main=paste(unique(as.character(tmp$Scaffold)),": ",round(diff(range(tmp$Consensus)),digit=1)," cM",sep=""),xlab="physical position (bp)",ylab="Genetic position (cM)",las=1,lty=1,lwd=2)
  ygene<-min(tmp[,3:6],na.rm=T)-5
  ygap<-c(min(tmp[,3:6],na.rm=T)-7,max(tmp[,3:6],na.rm=T))
  for(j in 1:dim(tmp1)[1]){
    if(tmp1$V3[j] == "gene"){
      lines(x=c(as.numeric(tmp1[j,3]),as.numeric(tmp1[j,4])),y=c(ygene,ygene),col="black",lwd=1,lty=1)
    }else{
      lines(x=c(as.numeric(tmp1[j,3]),as.numeric(tmp1[j,4])),y=c(ygene,ygene),col="black",lwd=3,lty=1)
    }
  }
  if(dim(tmp2)[1] > 0){
    for(k in 1:dim(tmp2)[1]){
      rect(tmp2$V2[k],ygap[1],tmp2$V3[k],ygap[2],density=NA,col=rgb(0.5,0.5,0.5,0.5))
    }
  } 
  points(tmp$Consensus_pos,tmp$Consensus,pch=19,lty=1,cex=1,col="black",type="b")
  points(tmp$C3_pos,tmp$C3,pch=17,lty=2,cex=1,col="green",type="b")
  points(tmp$C2_pos,tmp$C2,pch=18,lty=3,cex=1,col="red",type="b")
  points(tmp$C1_pos,tmp$C1,pch=6,lty=4,cex=1,col="blue",type="b")
}
dev.off()



multi_marker_scafs<-NULL
for (i in 1:length(levels(maps$Scaffold))){
  tmp <- maps[maps$Scaffold == levels(maps$Scaffold)[i],]
  if (dim(tmp)[1] > 1){
    multi_marker_scafs<-rbind(multi_marker_scafs,tmp)
  }
}

library(gdata)
multi_marker_scafs<-drop.levels(multi_marker_scafs)

split_scaf<-NULL
for(i in 1:length(levels(multi_marker_scafs$Scaffold))){
  split_scaf<-c(split_scaf,length(unique(multi_marker_scafs$LG[multi_marker_scafs$Scaffold == levels(multi_marker_scafs$Scaffold)[i]])))
}

pdf("Map_composition.pdf")
par(mfrow=c(1,1))
plot(table(table(maps$Scaffold))/length(levels(maps$Scaffold)),type="b",pch=19,ylim=c(0,1),las=1,xlab="Number of Markers",ylab="Fraction of Scaffolds")
library(Hmisc)
inner_plot<-subplot(plot(table(split_scaf)/length(levels(maps$Scaffold)),type="b",pch=19,las=1,ylab="",xlab="Number of LGs"),c(6,11),c(0.4,1))

op <- par(no.readonly=TRUE)
par(inner_plot)
points(1,dim(intra_split)[1]/length(levels(maps$Scaffold)),pch=19,col="red")
par(op)
dev.off()

rm(op,tmp,tmp2,i,j,k)


###Gene model confidance in split scaffolds
read.table("HC_gene-model_scaffolds.txt",head=F)->HC
read.table("MC_gene-model_scaffolds.txt",head=F)->MC
read.table("LC_gene-model_scaffolds.txt",head=F)->LC


### Which map scaffolds and how many scaffolds contain the different confidence level gene models 
HC_scaf_in_map<-maps$Scaffold[levels(maps$Scaffold) %in% levels(HC$V1)]
MC_scaf_in_map<-maps$Scaffold[levels(maps$Scaffold) %in% levels(MC$V1)]
LC_scaf_in_map<-maps$Scaffold[levels(maps$Scaffold) %in% levels(LC$V1)]

library(gdata)
HC_scaf_in_map<-droplevels(HC_scaf_in_map)
MC_scaf_in_map<-droplevels(MC_scaf_in_map)
LC_scaf_in_map<-droplevels(LC_scaf_in_map)


### Which confidence level of gene models are located on the split scaffolds?
#Inter split
length(inter_split$V1[inter_split$V1 %in% HC$V1 & inter_split$V1 %in% MC$V1 & inter_split$V1 %in% LC$V1])
length(inter_split$V1[inter_split$V1 %in% HC$V1 & inter_split$V1 %in% MC$V1 & !(inter_split$V1 %in% LC$V1)])
length(inter_split$V1[inter_split$V1 %in% HC$V1 & !(inter_split$V1 %in% MC$V1) & inter_split$V1 %in% LC$V1])
length(inter_split$V1[!(inter_split$V1 %in% HC$V1) & inter_split$V1 %in% MC$V1 & inter_split$V1 %in% LC$V1])
length(inter_split$V1[inter_split$V1 %in% HC$V1 & !(inter_split$V1 %in% MC$V1) & !(inter_split$V1 %in% LC$V1)])
length(inter_split$V1[!(inter_split$V1 %in% HC$V1) & inter_split$V1 %in% MC$V1 & !(inter_split$V1 %in% LC$V1)])
length(inter_split$V1[!(inter_split$V1 %in% HC$V1) & !(inter_split$V1 %in% MC$V1) & inter_split$V1 %in% LC$V1])
#Intra split
length(intra_split$V1[intra_split$V1 %in% HC$V1 & intra_split$V1 %in% MC$V1 & intra_split$V1 %in% LC$V1])
length(intra_split$V1[intra_split$V1 %in% HC$V1 & intra_split$V1 %in% MC$V1 & !(intra_split$V1 %in% LC$V1)])
length(intra_split$V1[intra_split$V1 %in% HC$V1 & !(intra_split$V1 %in% MC$V1) & intra_split$V1 %in% LC$V1])
length(intra_split$V1[!(intra_split$V1 %in% HC$V1) & intra_split$V1 %in% MC$V1 & intra_split$V1 %in% LC$V1])
length(intra_split$V1[intra_split$V1 %in% HC$V1 & !(intra_split$V1 %in% MC$V1) & !(intra_split$V1 %in% LC$V1)])
length(intra_split$V1[!(intra_split$V1 %in% HC$V1) & intra_split$V1 %in% MC$V1 & !(intra_split$V1 %in% LC$V1)])
length(intra_split$V1[!(intra_split$V1 %in% HC$V1) & !(intra_split$V1 %in% MC$V1) & intra_split$V1 %in% LC$V1])

save.image("split_scafs.RData")
