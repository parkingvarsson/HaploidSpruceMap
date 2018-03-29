#### Create Figure showing the number of markers per scaffold and how many of the scaffolds that are wrongly assembled
read.table("Consensus-maps.txt",head=T)->maps
maps$scaf<-as.factor(t(as.data.frame(strsplit(as.character(maps$marker),":")))[,1])

multi_marker_scafs<-NULL
for (i in 1:length(levels(maps$scaf))){
  tmp <- maps[maps$scaf == levels(maps$scaf)[i],]
  if (dim(tmp)[1] > 1){
    multi_marker_scafs<-rbind(multi_marker_scafs,tmp)
  }
}

library(gdata)
multi_marker_scafs<-drop.levels(multi_marker_scafs)

read.table("inter_split_scaffolds.txt",head=F)->inter_split
read.table("intra_split_scaffolds.txt",head=F)->intra_split


split_scaf<-NULL
for(i in 1:length(levels(multi_marker_scafs$scaf))){
  split_scaf<-c(split_scaf,length(unique(multi_marker_scafs$LG[multi_marker_scafs$scaf == levels(multi_marker_scafs$scaf)[i]])))
}

plot(table(table(maps$scaf))/length(levels(maps$scaf)),type="b",pch=19,ylim=c(0,1),las=1,xlab="Number of Markers",ylab="Fraction of Scaffolds")
library(Hmisc)
inner_plot<-subplot(plot(table(split_scaf)/length(levels(maps$scaf)),type="b",pch=19,las=1,ylab="",xlab="Number of LGs"),c(4,8),c(0.4,1))

op <- par(no.readonly=TRUE)
par(inner_plot)
points(1,dim(intra_split)[1]/length(levels(maps$scaf)),pch=19,col="red")
par(op)


###Gene model confidance in split scaffolds
read.table("HC_gene-model_scaffolds.txt",head=F)->HC
read.table("MC_gene-model_scaffolds.txt",head=F)->MC
read.table("LC_gene-model_scaffolds.txt",head=F)->LC


### Which map scaffolds and how many scaffolds contain the different confidence level gene models 
HC_scaf_in_map<-maps$scaf[levels(maps$scaf) %in% levels(HC$V1)]
MC_scaf_in_map<-maps$scaf[levels(maps$scaf) %in% levels(MC$V1)]
LC_scaf_in_map<-maps$scaf[levels(maps$scaf) %in% levels(LC$V1)]

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
