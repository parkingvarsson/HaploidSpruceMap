setwd("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/Cluster_comparison/dropped_maps/UPSCb_Consensus")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/Cluster_comparison/dropped_maps/UPSCb_Consensus/split_scafs.RData")

read.table("../Consensus_maps.txt",head=T)->maps

read.table("intra_split_scaffolds.txt",head=F)->intra_split
read.table("inter_split_scaffolds.txt",head=F)->inter_split

library(gdata)
inter_split_scaf<-drop.levels(maps[maps$Scaffold %in% inter_split$V1,])
#intra_split<-as.data.frame(drop.levels(intra_split[!(intra_split$V1 %in% c("MA_9458","MA_281725","MA_10431182","MA_10431315","MA_10432328")),]))
intra_split_scaf<-drop.levels(maps[maps$Scaffold %in% intra_split$V1,])

for (i in 1:length(levels(intra_split_scaf$Scaffold))){
  print(c(as.character(unique(intra_split_scaf$Scaffold[intra_split_scaf$Scaffold==levels(intra_split_scaf$Scaffold)[i]])),dim(intra_split_scaf[intra_split_scaf$Scaffold==levels(intra_split_scaf$Scaffold)[i],])[1],length(unique(intra_split_scaf$LG[intra_split_scaf$Scaffold==levels(intra_split_scaf$Scaffold)[i]]))))
}

for (i in 1:length(levels(inter_split_scaf$Scaffold))){
  print(c(as.character(unique(inter_split_scaf$Scaffold[inter_split_scaf$Scaffold==levels(inter_split_scaf$Scaffold)[i]])),dim(inter_split_scaf[inter_split_scaf$Scaffold==levels(inter_split_scaf$Scaffold)[i],])[1],length(unique(inter_split_scaf$LG[inter_split_scaf$Scaffold==levels(inter_split_scaf$Scaffold)[i]]))))
}

inter_split_scaf$Probe<-as.factor(inter_split_scaf$Probe)
inter_split_scaf$LG<-as.factor(inter_split_scaf$LG)

inter_split_structure<-data.frame(matrix(ncol = 13, nrow = 164))
names(inter_split_structure)<-c("Scaffold","LG I","LG II","LG III","LG IV","LG V","LG VI","LG VII","LG VIII","LG IX","LG X","LG XI","LG XII")
for(i in 1:length(levels(inter_split_scaf$Scaffold))){
  inter_split_structure[i,1]<-levels(inter_split_scaf$Scaffold)[i]
  for(j in 1:length(levels(inter_split_scaf$LG))){
    if (length(unique(inter_split_scaf$Probe[inter_split_scaf$Scaffold == levels(inter_split_scaf$Scaffold)[i] & inter_split_scaf$LG==levels(inter_split_scaf$LG)[j]])) >0){
      tmp<-as.character(unique(inter_split_scaf$Probe[inter_split_scaf$Scaffold == levels(inter_split_scaf$Scaffold)[i] & inter_split_scaf$LG==levels(inter_split_scaf$LG)[j]]))
      if(length(tmp)==1){
        inter_split_structure[i,j+1]<-tmp
      }
      if(length(tmp)==2){
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],sep=",")
      }
      if(length(tmp)==3){
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],tmp[3],sep=",")
      }
      if(length(tmp)==4){
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],tmp[3],tmp[4],sep=",")
      }
      if(length(tmp)==5){
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],sep=",")
      }
      if(length(tmp)==6){
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],sep=",")
      }
      if(length(tmp)==7){
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],sep=",")
      }
      if(length(tmp)==8){
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8],sep=",")
      }
      if(length(tmp)==9){
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8],tmp[9],sep=",")
      }
      if(length(tmp)==10){
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8],tmp[9],tmp[10],sep=",")
      }
    }
    else{ 
      inter_split_structure[i,j+1]<-"-"
    }  
  }
}

write.table(inter_split_structure,"Structure_of_split_scaffolds.txt",sep="\t",quote=F,row.names =F)



extended_probes<-read.table("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/extended_probe_positions.bed.txt",sep="\t",head=T)
extended_probes<-extended_probes[do.call(order,extended_probes),]
scaf<-NULL
LG<-NULL
start<-NULL
end<-NULL
for (i in 1:dim(inter_split_structure)[1]){
  for (j in 1:12){
    if (inter_split_structure[i,j+1] != "-"){
      df<-as.numeric(strsplit(inter_split_structure[i,j+1],",")[[1]])
      for(k in 1:length(df)){
        scaf<-c(scaf,inter_split_structure$Scaffold[i])
        LG<- c(LG,j)
        start<-c(start,extended_probes[extended_probes$CHROM == inter_split_structure$Scaffold[i],][df[k],2])
        end<-c(end,extended_probes[extended_probes$CHROM == inter_split_structure$Scaffold[i],][df[k],3])
      }
    }  
  }
}
rm(df)
split_position<-as.data.frame(cbind(scaf,LG,start,end))
write.table(split_position,"split_position.txt",quote=F,row.names = F,sep="\t")


## probe positions in scaffolds that are split between LGs
messed_up_scafs<-NULL
for (i in 1:length(levels(split_position$scaf))){
  tmp<-split_position[split_position$scaf== levels(split_position$scaf)[i],]
  tmp$start<-as.integer(as.character(tmp$start))
  tmp$end<-as.integer(as.character(tmp$end))
  for (j in 1:(dim(tmp)[1]-1)){
    for (m in 2:dim(tmp)[1]){
      if ((tmp[m,3] < tmp[j,3]  & tmp[j,3] < tmp[m,4]) | (tmp[m,3]< tmp[j,4] & tmp[j,4] < tmp[m,4]) | (tmp[j,3] < tmp[m,3]  & tmp[m,3] < tmp[j,4]) | (tmp[j,3]< tmp[m,4] & tmp[m,4] < tmp[j,4])){
        messed_up_scafs<-c(messed_up_scafs,as.character(tmp$scaf[1]))
      }
    }
  }
}


split_probes<-extended_probes[extended_probes$CHROM %in% messed_up_scafs,]
library(gdata)
split_probes<-drop.levels(split_probes)

read.table("haploid_scaffolds.txt",head=F)->scaf_length
read.table("HC_genes.gff3",head=F)->HC
read.table("MC_genes.gff3",head=F)->MC
read.table("LC_genes.gff3",head=F)->LC
read.table("Pabies1.0-all.phase.changed.gff3",head=F)->All_models
All_models<-drop.levels(All_models[All_models$V1 %in% as.character(inter_split$V1),])
read.table("inter_split_Scaffolds.gaps.bed",head=F)->gaps

pdf("Split_positions-InterSplit_scaffolds.pdf")
par(mfrow=c(3,3))
for(i in 1:length(levels(split_position$scaf))){
  tmp<-split_position[split_position$scaf==levels(split_position$scaf)[i],]
  tmp1<-HC[as.character(HC$V1)==as.character(unique(tmp$scaf)),c(1,4,5)]
  tmp2<-MC[as.character(MC$V1)==as.character(unique(tmp$scaf)),c(1,4,5)]
  tmp3<-LC[as.character(LC$V1)==as.character(unique(tmp$scaf)),c(1,4,5)]
  tmp4<-All_models[as.character(All_models$V1)==as.character(unique(tmp$scaf)) & All_models$V3 %in% c("gene","exon"),c(1,3,4,5)]
  tmp5<-gaps[as.character(gaps$V1)==as.character(unique(tmp$scaf)),]
  
  plot(1,type="n",ylim=c(-0.5,12.5),xlim=c(0,scaf_length$V2[scaf_length$V1==as.character(unique(tmp$scaf))]),las=1,xlab="Scaffold length (bp)",ylab="LG",main=as.character(unique(tmp$scaf)),axes=F)
  abline(h=c(seq(1,12,1)),col="gray80",lwd=0.5)
  axis(1);box()
  mtext(side=2,at=c(seq(0,12,1)),text=c("Gene",seq(1,12,1)),las=1,cex=0.5,line=1)
  for(j in 1:dim(tmp)[1]){
    lines(x=c(as.numeric(as.character(tmp$start[j])),as.numeric(as.character(tmp$end[j]))),y=c(as.numeric(as.character(tmp$LG[j])),as.numeric(as.character(tmp$LG[j]))),col="red",lwd=3)
  }
#  if(dim(tmp1)[1] > 0){
#    for(k in 1:dim(tmp1)[1]){
#    lines(x=c(as.numeric(tmp1[k,2]),as.numeric(tmp1[k,3])),y=c(0,0),col="black",lwd=1,lty=1)
#    }
#  }  
#  if(dim(tmp2)[1] > 0){
#    for(l in 1:dim(tmp2)[1]){
#    lines(x=c(as.numeric(tmp2[l,2]),as.numeric(tmp2[l,3])),y=c(-1,-1),col="black",lwd=1,lty=1)
#    }
#  }  
#  if(dim(tmp3)[1] > 0){
#    for(m in 1:dim(tmp3)[1]){
#    lines(x=c(as.numeric(tmp3[m,2]),as.numeric(tmp3[m,3])),y=c(-2,-2),col="black",lwd=1,lty=1)
#    }
#  }
  for(n in 1:dim(tmp4)[1]){
    if(tmp4$V3[n] == "gene"){
      lines(x=c(as.numeric(tmp4[n,3]),as.numeric(tmp4[n,4])),y=c(-0,-0),col="black",lwd=1,lty=1)
    }else{
      lines(x=c(as.numeric(tmp4[n,3]),as.numeric(tmp4[n,4])),y=c(-0,-0),col="black",lwd=3,lty=1)
    }
  }
  if(dim(tmp5)[1] > 0){
    for(o in 1:dim(tmp5)[1]){
      rect(tmp5$V2[o],-0.5,tmp5$V3[o],12.5,density=NA,col=rgb(0.5,0.5,0.5,0.5))
    }
  }  
}
dev.off()
rm(tmp,tmp1,tmp2,tmp3,tmp4,tmp5)

save.image("Split_Scaffold-structure.RData")

multi_marker_scafs<-NULL
for (i in 1:length(levels(maps$Scaffold))){
  tmp <- maps[maps$Scaffold == levels(maps$Scaffold)[i],]
  if (dim(tmp)[1] > 1){
    multi_marker_scafs<-rbind(multi_marker_scafs,tmp)
  }
}

multi_marker_scafs<-drop.levels(multi_marker_scafs)




