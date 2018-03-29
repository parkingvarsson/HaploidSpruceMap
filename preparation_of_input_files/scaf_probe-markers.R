##probe info
read.table("extended_probe_positions.txt",head=T)->probes
probes[do.call(order,probes),]->probes
probes$probe<-NULL
probes$scaf_probe<-NULL
for (i in 1:length(levels(probes$CHROM))){
  for (j in 1:dim(probes[probes$CHROM == levels(probes$CHROM)[i],])[1]){
    probes$probe[probes$CHROM == levels(probes$CHROM)[i]][j]<-paste("p",j,sep="")
    probes$scaf_probe[probes$CHROM == levels(probes$CHROM)[i]][j]<-paste(levels(probes$CHROM)[i],j,sep=":")
  }
}


#set which scaf_probe each marker belong to in the different clusters
findProbe<-function(probes,info_markers){
  df<-probes$scaf_probe[probes$CHROM %in% info_markers$CHROM & probes$Start-1 < info_markers$POS & probes$End+1 > info_markers$POS]
  return(df)
}


info_markers_C1$scaf_probe<-NULL
for (i in 1:dim(info_markers_C1)[1]){
  info_markers_C1$scaf_probe[i]<-findProbe(probes,info_markers_C1[i,])
}
info_markers_C1$scaf_probe<-as.factor(info_markers_C1$scaf_probe)

info_markers_C2$scaf_probe<-NULL
for (i in 1:dim(info_markers_C2)[1]){
  info_markers_C2$scaf_probe[i]<-findProbe(probes,info_markers_C2[i,])
}
info_markers_C2$scaf_probe<-as.factor(info_markers_C2$scaf_probe)

info_markers_C3$scaf_probe<-NULL
for (i in 1:dim(info_markers_C3)[1]){
  info_markers_C3$scaf_probe[i]<-findProbe(probes,info_markers_C3[i,])
}
info_markers_C3$scaf_probe<-as.factor(info_markers_C3$scaf_probe)

## create a frequency table of how many clusters each scaf_probe exists in
scaf_probe_occurance<-as.data.frame(table(c(levels(info_markers_C1$scaf_probe),levels(info_markers_C2$scaf_probe),levels(info_markers_C3$scaf_probe))))

venn(list(levels(info_markers_C1$scaf_probe),levels(info_markers_C2$scaf_probe),levels(info_markers_C3$scaf_probe)))