### Read in which scaffolds that contain HC, MC and LC gene models
read.table("HC_gene-model_scaffolds.txt",head=F)->HC
read.table("MC_gene-model_scaffolds.txt",head=F)->MC
read.table("LC_gene-model_scaffolds.txt",head=F)->LC

### Read in the consensus haploid map 
read.table("Consensus-maps.txt",head=T)->map

## Extract the scaffold name from the marker names and put them in an extra column
map$scaf<-NULL
for (i in 1:dim(map)[1]){
  map$scaf[i]<-strsplit(as.character(map$marker[i]),split=":")[[1]][1]
}
map$scaf<-as.factor(map$scaf)

### Which map scaffolds and how many scaffolds contain the different confidence level gene models 
HC_scaf_in_map<-map$scaf[levels(map$scaf) %in% levels(HC$V1)]
MC_scaf_in_map<-map$scaf[levels(map$scaf) %in% levels(MC$V1)]
LC_scaf_in_map<-map$scaf[levels(map$scaf) %in% levels(LC$V1)]

library(gdata)
HC_scaf_in_map<-droplevels(HC_scaf_in_map)
MC_scaf_in_map<-droplevels(MC_scaf_in_map)
LC_scaf_in_map<-droplevels(LC_scaf_in_map)


### Read in the inter and intra split scaffold information
read.table("intra_split_scaffolds.txt",head=F)->intra_split
read.table("inter_split_scaffolds.txt",head=F)->inter_split

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





