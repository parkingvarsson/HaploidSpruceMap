##Calculate Zns (Kelly 1997) per probe using bed file to determine probe start and stop

library("data.table")
fread("r2_320bp_windows.geno.ld",head=T,stringsAsFactors=F) -> ldout
names(ldout)[5]<-"R2"
fread("../overlapping_probes_with_map_new_with_probe_and_scaffold_info.bed",head=T,stringsAsFactors=F) -> bed

Z<-NULL
for (i in 1:dim(bed)[1]){
    df<-subset(ldout,CHR==bed$scaffold[i] & POS1 > bed$start[i] & POS1 < bed$stop[i] & POS2 > bed$start[i] & POS2 < bed$stop[i])
  
  tmp<-NULL
  tmp$scaffold<-bed$scaffold[i]
  tmp$probe<-bed$scaf_prob[i]
  tmp$start<-min(c(df$POS1,df$POS2))
  tmp$stop<-max(c(df$POS1,df$POS2))
  tmp$mid<-mean(c(tmp$start,tmp$stop))
  tmp$SNPs<-length(unique(c(df$POS1,df$POS2)))
  tmp$r2<-sum(df$R2,na.rm=T)
  tmp$z<-(2*tmp$r2)/(tmp$SNPs*(tmp$SNPs-1))
  tmp<-as.data.frame(tmp)
  Z<-rbind(Z,tmp)
}
Z<-as.data.frame(Z)
write.table(Z,"Zns_per_probe.txt",quote=F,row.names=F,sep="\t")
