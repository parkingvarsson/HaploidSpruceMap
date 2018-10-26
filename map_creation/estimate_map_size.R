load("/Users/caabon02/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_1/C1_cleaned_maps.RData")

library(BatchMap)

estimate_size<-NULL
estimate_size$LG<-NULL
estimate_size$size<-NULL
for(i in 1:12){
  for(j in 1:100){
    print(c(i,j))
    df<-map(record(make.seq(twopt_clean,sample(LGs[[i]][[1]],100,replace = F))),set.map.fun("kosambi")) 
    estimate_size$LG<-c(estimate_size$LG,i)
    estimate_size$size<- c(estimate_size$size,max(cumsum(c(0, get(get(".map.fun", envir = .onemapEnv))(df$seq.rf)))))
  }
}  
estimate_size<-as.data.frame(estimate_size)

write.table(estimate_size,"C1_estimated_size.txt",row.names = F,col.names = T,quote=F, sep="\t")
save(estimate_size,file="estimate_size.RData")