read.table("Consensus-maps.txt",head=T)->maps

library(gdata)
inter_split_scaf<-drop.levels(maps[maps$scaf %in% inter_split$V1,])
intra_split_scaf<-drop.levels(maps[maps$scaf %in% intra_split$V1,])

for (i in 1:length(levels(intra_split_scaf$scaf))){
  print(c(as.character(unique(intra_split_scaf$scaf[intra_split_scaf$scaf==levels(intra_split_scaf$scaf)[i]])),dim(intra_split_scaf[intra_split_scaf$scaf==levels(intra_split_scaf$scaf)[i],])[1],length(unique(intra_split_scaf$LG[intra_split_scaf$scaf==levels(intra_split_scaf$scaf)[i]]))))
}

for (i in 1:length(levels(inter_split_scaf$scaf))){
  print(c(as.character(unique(inter_split_scaf$scaf[inter_split_scaf$scaf==levels(inter_split_scaf$scaf)[i]])),dim(inter_split_scaf[inter_split_scaf$scaf==levels(inter_split_scaf$scaf)[i],])[1],length(unique(inter_split_scaf$LG[inter_split_scaf$scaf==levels(inter_split_scaf$scaf)[i]]))))
}

inter_split_scaf$probe<-as.factor(t(as.data.frame(strsplit(as.character(inter_split_scaf$marker),":")))[,2])
inter_split_scaf$LG<-as.factor(inter_split_scaf$LG)

inter_split_structure<-data.frame(matrix(ncol = 13, nrow = 164))
names(inter_split_structure)<-c("Scaffold","LG I","LG II","LG III","LG IV","LG V","LG VI","LG VII","LG VIII","LG IX","LG X","LG XI","LG XII")
for(i in 1:length(levels(inter_split_scaf$scaf))){
  inter_split_structure[i,1]<-levels(inter_split_scaf$scaf)[i]
  for(j in 1:length(levels(inter_split_scaf$LG))){
    if (length(unique(inter_split_scaf$probe[inter_split_scaf$scaf == levels(inter_split_scaf$scaf)[i] & inter_split_scaf$LG==levels(inter_split_scaf$LG)[j]])) >0){
      tmp<-as.character(unique(inter_split_scaf$probe[inter_split_scaf$scaf == levels(inter_split_scaf$scaf)[i] & inter_split_scaf$LG==levels(inter_split_scaf$LG)[j]]))
      if(length(tmp)==1){
        inter_split_structure[i,j+1]<-tmp
      }else{
        inter_split_structure[i,j+1]<-paste(tmp[1],tmp[2],sep=",")
      }
    }
    else{ 
      inter_split_structure[i,j+1]<-"-"
    }  
  }
}


## Read in table with spplit scaf before marker dropping
split_before_drop<-read.table("Splitted_scaffolds.txt",sep="\t",stringsAsFactor=F)
split_before_drop<-split_before_drop[,c(1:5,7:10,6,11:13)]
names(split_before_drop)<-c("Scaffold","LG I","LG II","LG III","LG IV","LG V","LG VI","LG VII","LG VIII","LG IX","LG X","LG XI","LG XII")
write.table(split_before_drop,"Split_scaf_pre-drop_right-names.txt",sep="\t",quote=F,row.names =F)

full_inter_split_structure<-rbind(split_before_drop,inter_split_structure[!(inter_split_structure$Scaffold %in% split_before_drop$Scaffold),])
full_inter_split_structure<-full_inter_split_structure[order(full_inter_split_structure$Scaffold),]
write.table(full_inter_split_structure,"Structure_of_split_scaffolds.txt",sep="\t",quote=F,row.names =F)

extended_probes<-read.table("extended_probe_positions.txt",sep="\t",head=T)
extended_probes<-extended_probes[do.call(order,extended_probes),]
scaf<-NULL
LG<-NULL
start<-NULL
end<-NULL
for (i in 1:dim(full_inter_split_structure)[1]){
  for (j in 1:12){
    if (full_inter_split_structure[i,j+1] != "-"){
      scaf<-c(scaf,full_inter_split_structure$Scaffold[i])
      LG<- c(LG,j)
      start<-c(start,extended_probes[extended_probes$CHROM == full_inter_split_structure$Scaffold[i],][min(as.numeric(strsplit(full_inter_split_structure[i,j+1],",")[[1]])),2])
      end<-c(end,extended_probes[extended_probes$CHROM == full_inter_split_structure$Scaffold[i],][max(as.numeric(strsplit(full_inter_split_structure[i,j+1],",")[[1]])),3])
    }
  }
}

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









