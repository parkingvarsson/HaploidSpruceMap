read.table("Consensus-maps.txt",head=T)->maps
read.table("HC_genes.gff3",head=F)->HC
read.table("MC_genes.gff3",head=F)->MC
read.table("LC_genes.gff3",head=F)->LC
names(HC)<-c("scaffold","Master","gene","start","end","not_known","strand","not_known2","gene_model")
names(MC)<-c("scaffold","Master","gene","start","end","not_known","strand","not_known2","gene_model")
names(LC)<-c("scaffold","Master","gene","start","end","not_known","strand","not_known2","gene_model")

maps$scaf<-as.factor(t(as.data.frame(strsplit(as.character(maps$marker),":")))[,1])
HC$gene_id<-as.factor(t(as.data.frame(strsplit(as.character(HC$gene_model),"=")))[,3])
MC$gene_id<-as.factor(t(as.data.frame(strsplit(as.character(MC$gene_model),"=")))[,3])
LC$gene_id<-as.factor(t(as.data.frame(strsplit(as.character(LC$gene_model),"=")))[,3])

###Number of gene model in mapped scaffolds
sum(c(dim(HC[HC$scaffold %in% maps$scaf,])[1],dim(MC[MC$scaffold %in% maps$scaf,])[1],dim(LC[LC$scaffold %in% maps$scaf,])[1]))
#17079
###Number of gene models in total
sum(c(dim(HC)[1],dim(MC)[1],dim(LC)[1]))
#66632

### % of gene models in mapped scaffolds
17079/66632
#25.63%

###HC models
dim(HC[HC$scaffold %in% maps$scaf,])[1]/dim(HC)[1]
# 31.69% (8379 out of 26437 in mapped scaffolds)

###MC models
dim(MC[MC$scaffold %in% maps$scaf,])[1]/dim(MC)[1]
# 20.6% (6624 out of 32150 in mapped scaffolds)

###MC models
dim(LC[LC$scaffold %in% maps$scaf,])[1]/dim(LC)[1]
# 25.8% (2076 out of 8045 in mapped scaffolds)


read.table("inter_split_scaffolds.txt",head=F)->inter_split_scaf
read.table("intra_split_scaffolds.txt",head=F)->intra_split_scaf

dim(HC[HC$scaffold %in% intra_split_scaf$V1,])[1] #15
dim(MC[MC$scaffold %in% intra_split_scaf$V1,])[1] #15
dim(LC[LC$scaffold %in% intra_split_scaf$V1,])[1] #6

dim(HC[HC$scaffold %in% inter_split_scaf$V1,])[1] #145
dim(MC[MC$scaffold %in% inter_split_scaf$V1,])[1] #114
dim(LC[LC$scaffold %in% inter_split_scaf$V1,])[1] #35
###330 gene models in total are located on the split scaffolds

dim(HC[HC$scaffold %in% as.factor(names(table(maps$scaf)[table(maps$scaf)>1])),])[1]
dim(MC[MC$scaffold %in% as.factor(names(table(maps$scaf)[table(maps$scaf)>1])),])[1]
dim(LC[LC$scaffold %in% as.factor(names(table(maps$scaf)[table(maps$scaf)>1])),])[1]

###Number of gene models on multi-marker scaffolds before marker were dropped
read.table("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/Cluster_comparison/Cluster1-all_rec40_maps.txt",head=F)->C1
read.table("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/Cluster_comparison/Cluster2-all_rec40_maps.txt",head=F)->C2
read.table("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/Cluster_comparison/Cluster3_40rec-Maps.txt",head=F)->C3

unique(c(as.character(C1$V2),as.character(C2$V2),as.character(C3$V2)))->unique_scafProbes
as.factor(t(as.data.frame(strsplit(unique_scafProbes,":")))[,1])->Scaffolds

dim(HC[HC$scaffold %in% as.factor(names(table(Scaffolds)[table(Scaffolds)>1])),])[1]
dim(MC[MC$scaffold %in% as.factor(names(table(Scaffolds)[table(Scaffolds)>1])),])[1]
dim(LC[LC$scaffold %in% as.factor(names(table(Scaffolds)[table(Scaffolds)>1])),])[1]







