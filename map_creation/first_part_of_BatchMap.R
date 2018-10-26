setwd("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_1")

library(BatchMap)
C1.out<-read.outcross2("Outcross_C1_probe_new-version.txt")

bins<-find.bins(C1.out,exact=F)
C1.out_clean<-create.data.bins(C1.out,bins)

twopt_clean<-rf.2pts(C1.out_clean,LOD=8,max.rf=0.35)
LGs_clean<-group(make.seq(twopt_clean,"all"))
print(LGs_clean,detailed=F)

set.map.fun(type="kosambi")

for(i in 1:12){
  print(i)
  assign(paste("C1_LG",i,sep="_"),make.seq(LGs_clean,i))
}

LGs<-list(C1_LG_1=C1_LG_1,C1_LG_2=C1_LG_2,C1_LG_3=C1_LG_3,C1_LG_4=C1_LG_4,C1_LG_5=C1_LG_5,C1_LG_6=C1_LG_6,C1_LG_7=C1_LG_7,C1_LG_8=C1_LG_8,C1_LG_9=C1_LG_9,C1_LG_10=C1_LG_10,C1_LG_11=C1_LG_11,C1_LG_12=C1_LG_12)

rm(C1_LG_1,C1_LG_2,C1_LG_3,C1_LG_4,C1_LG_5,C1_LG_6,C1_LG_7,C1_LG_8,C1_LG_9,C1_LG_10,C1_LG_11,C1_LG_12,i)

save.image("C1_cleaned.RData")

