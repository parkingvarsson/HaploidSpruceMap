
Cluster1<-read.table("Cluster_1-no-cont.GT.FORMAT.txt",head=T)
Cluster2<-read.table("Cluster_2-no-cont.GT.FORMAT.txt",head=T)
Cluster3<-read.table("Cluster_3-no-cont.GT.FORMAT.txt",head=T)

### covert strange genotypes into missing data (only keep ./., 1/1, 0/0)
for (i in 3:dim(Cluster1)[2]){
  Cluster1[,i]<-as.factor(gsub("2/2","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("0/2","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("0/3","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("1/2","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("2/1","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("3/3","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("3/2","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("1/3","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("2/3","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("0/1","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("4/4","./.",Cluster1[,i]))
  Cluster1[,i]<-as.factor(gsub("0/4","./.",Cluster1[,i]))
}


for (i in 3:dim(Cluster2)[2]){
  Cluster2[,i]<-as.factor(gsub("2/2","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("0/2","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("0/3","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("1/2","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("2/1","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("3/3","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("3/2","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("1/3","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("2/3","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("0/1","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("4/4","./.",Cluster2[,i]))
  Cluster2[,i]<-as.factor(gsub("0/4","./.",Cluster2[,i]))
}


for (i in 3:dim(Cluster3)[2]){
  Cluster3[,i]<-as.factor(gsub("2/2","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("0/2","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("0/3","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("1/2","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("2/1","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("3/3","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("3/2","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("1/3","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("2/3","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("0/1","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("4/4","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("0/4","./.",Cluster3[,i]))
  Cluster3[,i]<-as.factor(gsub("3/1","./.",Cluster3[,i]))
}



## Subset markers with less than 20% missing data (Cluster1 < 62; Cluster2 < 54; Cluster 3 < 168) 
tmp<-apply(Cluster1[,3:316],1,function(x){length(table(x))>2})
Cluster1_red<-Cluster1[which(tmp=="TRUE"),]
Cluster1_red_table<-as.data.frame(apply(Cluster1_red[,3:316],1,table))
Cluster1_final<-Cluster1_red[which(Cluster1_red_table[1,]<62),]
Cluster1_final<-drop.levels(Cluster1_final)

tmp<-apply(Cluster2[,3:272],1,function(x){length(table(x))>2})
Cluster2_red<-Cluster2[which(tmp=="TRUE"),]
Cluster2_red_table<-as.data.frame(apply(Cluster2_red[,3:272],1,table))
Cluster2_final<-Cluster2_red[which(Cluster2_red_table[1,]<54),] 
Cluster2_final<-drop.levels(Cluster2_final)

tmp<-apply(Cluster3[,3:844],1,function(x){length(table(x))>2})
Cluster3_red<-Cluster3[which(tmp=="TRUE"),]
Cluster3_red_table<-as.data.frame(apply(Cluster3_red[,3:844],1,table))
Cluster3_final<-Cluster3_red[which(Cluster3_red_table[1,]<168),]
Cluster3_final<-drop.levels(Cluster3_final)



### recode missing data to "-"  and 0/0 -> "a" and 1/1 -> "b"
for (i in 3:dim(Cluster1_final)[2]){
  Cluster1_final[,i]<-gsub("0/0","a",Cluster1_final[,i])
  Cluster1_final[,i]<-gsub("1/1","b",Cluster1_final[,i])
  Cluster1_final[,i]<-gsub("./.","-",Cluster1_final[,i])
}

for (i in 3:dim(Cluster2_final)[2]){
  Cluster2_final[,i]<-gsub("0/0","a",Cluster2_final[,i])
  Cluster2_final[,i]<-gsub("1/1","b",Cluster2_final[,i])
  Cluster2_final[,i]<-gsub("./.","-",Cluster2_final[,i])
}

for (i in 3:dim(Cluster3_final)[2]){
  Cluster3_final[,i]<-gsub("0/0","a",Cluster3_final[,i])
  Cluster3_final[,i]<-gsub("1/1","b",Cluster3_final[,i])
  Cluster3_final[,i]<-gsub("./.","-",Cluster3_final[,i])
}

### merge GT data into 1 column and count alleles 
C1<-NULL
C1$CHROM<-Cluster1_final$CHROM
C1$POS<-Cluster1_final$POS
C1$GT<-NULL
for (i in 1:dim(Cluster1_final)[1]){
  C1$GT[i]<-paste(Cluster1_final[i,3:316],sep="",collapse=",")
}
C1<-as.data.frame(C1)

library(stringi)
for (i in 1:dim(C1)[1]){
  C1$miss[i]<-stri_count_fixed(C1$GT[i],c("a","b","-")[3])
  C1$a[i]<-stri_count_fixed(C1$GT[i],c("a","b","-")[1])
  C1$b[i]<-stri_count_fixed(C1$GT[i],c("a","b","-")[2])
}


C2<-NULL
C2$CHROM<-Cluster2_final$CHROM
C2$POS<-Cluster2_final$POS
C2$GT<-NULL
for (i in 1:dim(Cluster2_final)[1]){
  C2$GT[i]<-paste(Cluster2_final[i,3:272],sep="",collapse=",")
}
C2<-as.data.frame(C2)

library(stringi)
for (i in 1:dim(C2)[1]){
  C2$miss[i]<-stri_count_fixed(C2$GT[i],c("a","b","-")[3])
  C2$a[i]<-stri_count_fixed(C2$GT[i],c("a","b","-")[1])
  C2$b[i]<-stri_count_fixed(C2$GT[i],c("a","b","-")[2])
}

C3<-NULL
C3$CHROM<-Cluster3_final$CHROM
C3$POS<-Cluster3_final$POS
C3$GT<-NULL
for (i in 1:dim(Cluster3_final)[1]){
  C3$GT[i]<-paste(Cluster3_final[i,3:844],sep="",collapse=",")
}
C3<-as.data.frame(C3)

library(stringi)
for (i in 1:dim(C3)[1]){
  C3$miss[i]<-stri_count_fixed(C3$GT[i],c("a","b","-")[3])
  C3$a[i]<-stri_count_fixed(C3$GT[i],c("a","b","-")[1])
  C3$b[i]<-stri_count_fixed(C3$GT[i],c("a","b","-")[2])
}

### Function to calculate the minor allele frequence (MAF)
maf<-function(a,b){
  if(a>b){
    maf<-b/(a+b)
  }
  else{
    maf<-a/(a+b)
  }
  return (maf)  
}

C1$maf<-NULL
for (i in 1:dim(C1)[1]){
  C1$maf[i]<-maf(C1$a[i],C1$b[i])
}

C2$maf<-NULL
for (i in 1:dim(C2)[1]){
  C2$maf[i]<-maf(C2$a[i],C2$b[i])
}

C3$maf<-NULL
for (i in 1:dim(C3)[1]){
  C3$maf[i]<-maf(C3$a[i],C3$b[i])
}


### extract the informative markers per Cluster with maf >0.4 ### don't care about godness of fit test
info_markers_C1<-C1[C1$maf>0.4,]
info_markers_C2<-C2[C2$maf>0.4,]
info_markers_C3<-C3[C3$maf>0.4,]

info_markers_C1<-drop.levels(info_markers_C1)
info_markers_C2<-drop.levels(info_markers_C2)
info_markers_C3<-drop.levels(info_markers_C3)


to_onemap_C1<-NULL  
for (i in 1:length(unique(info_markers_C1$CHROM))){
  SNPs<-NULL
  if (dim(info_markers_C1[info_markers_C1$CHROM==unique(info_markers_C1$CHROM)[i],])[1]==1){
    to_onemap_C1<-rbind(to_onemap_C1,info_markers_C1[info_markers_C1$CHROM==unique(info_markers_C1$CHROM)[i],])
  }
  else {
    SNPs<- info_markers_C1[info_markers_C1$CHROM==unique(info_markers_C1$CHROM)[i],]
    SNPs<-SNPs[SNPs$miss==min(SNPs$miss),]
    SNPs<-SNPs[SNPs$maf==max(SNPs$maf),][1,]
    to_onemap_C1<-rbind(to_onemap_C1, SNPs)
  }
}

to_onemap_C2<-NULL  
for (i in 1:length(unique(info_markers_C2$CHROM))){
  SNPs<-NULL
  if (dim(info_markers_C2[info_markers_C2$CHROM==unique(info_markers_C2$CHROM)[i],])[1]==1){
    to_onemap_C2<-rbind(to_onemap_C2,info_markers_C2[info_markers_C2$CHROM==unique(info_markers_C2$CHROM)[i],])
  }
  else {
    SNPs<- info_markers_C2[info_markers_C2$CHROM==unique(info_markers_C2$CHROM)[i],]
    SNPs<-SNPs[SNPs$maf==max(SNPs$maf),][1,]
    to_onemap_C2<-rbind(to_onemap_C2, SNPs)
  }
}

to_onemap_C3<-NULL  
for (i in 1:length(unique(info_markers_C3$CHROM))){
  SNPs<-NULL
  if (dim(info_markers_C3[info_markers_C3$CHROM==unique(info_markers_C3$CHROM)[i],])[1]==1){
    to_onemap_C3<-rbind(to_onemap_C3,info_markers_C3[info_markers_C3$CHROM==unique(info_markers_C3$CHROM)[i],])
  }
  else {
    SNPs<- info_markers_C3[info_markers_C3$CHROM==unique(info_markers_C3$CHROM)[i],]
    SNPs<-SNPs[SNPs$maf==max(SNPs$maf),][1,]
    to_onemap_C3<-rbind(to_onemap_C3, SNPs)
  }
}



### convert trustworthy SNPs (from "to_onemap") to Onemap input format
Onemap_C1<-NULL
Onemap_C1$cross<-NULL
Onemap_C1$marker<-NULL
for (i in 1:dim(to_onemap_C1)[1]){
  Onemap_C1$marker[i]<-paste("*",to_onemap_C1$CHROM[i],sep="")
  Onemap_C1$cross[i]<-"D1.11"
}
Onemap_C1$GT<-to_onemap_C1$GT
Onemap_C1<-as.data.frame(Onemap_C1)


Onemap_C2<-NULL
Onemap_C2$cross<-NULL
Onemap_C2$marker<-NULL
for (i in 1:dim(to_onemap_C2)[1]){
  Onemap_C2$marker[i]<-paste("*",to_onemap_C2$CHROM[i],sep="",collapse="")
  Onemap_C2$cross[i]<-"D1.11"
}
Onemap_C2$GT<-to_onemap_C2$GT
Onemap_C2<-as.data.frame(Onemap_C2)


Onemap_C3<-NULL
Onemap_C3$cross<-NULL
Onemap_C3$marker<-NULL
for (i in 1:dim(to_onemap_C3)[1]){
  Onemap_C3$marker[i]<-paste("*",to_onemap_C3$CHROM[i],sep="")
  Onemap_C3$cross[i]<-"D1.11"
}
Onemap_C3$GT<-to_onemap_C3$GT
Onemap_C3<-as.data.frame(Onemap_C3)

save.image(".RData")

write.table(Onemap_C1,"known_SNPs_Onemap_C1.txt",row.names=F,col.names=F,quote=F)
write.table(Onemap_C2,"known_SNPs_Onemap_C2.txt",row.names=F,col.names=F,quote=F)
write.table(Onemap_C3,"known_SNPs_Onemap_C3.txt",row.names=F,col.names=F,quote=F)



## create Onemap input file for the clusters at probe level
##C1
to_onemap_C1_probe<-NULL  
for (i in 1:length(unique(info_markers_C1$scaf_probe))){
  SNPs<-NULL
  if (dim(info_markers_C1[info_markers_C1$scaf_probe==unique(info_markers_C1$scaf_probe)[i],])[1]==1){
    to_onemap_C1_probe<-rbind(to_onemap_C1_probe,info_markers_C1[info_markers_C1$scaf_probe==unique(info_markers_C1$scaf_probe)[i],])
  }
  else {
    SNPs<- info_markers_C1[info_markers_C1$scaf_probe==unique(info_markers_C1$scaf_probe)[i],]
    SNPs<-SNPs[SNPs$miss==min(SNPs$miss),]
    SNPs<-SNPs[SNPs$maf==max(SNPs$maf),][1,]
    to_onemap_C1_probe<-rbind(to_onemap_C1_probe, SNPs)
  }
}

Onemap_C1_probe<-NULL
Onemap_C1_probe$cross<-NULL
Onemap_C1_probe$marker<-NULL
for (i in 1:dim(to_onemap_C1_probe)[1]){
  Onemap_C1_probe$marker[i]<-paste("*",to_onemap_C1_probe$scaf_probe[i],sep="")
  Onemap_C1_probe$cross[i]<-"D1.11"
}
Onemap_C1_probe$GT<-to_onemap_C1_probe$GT
Onemap_C1_probe<-as.data.frame(Onemap_C1_probe)


##C2
to_onemap_C2_probe<-NULL  
for (i in 1:length(unique(info_markers_C2$scaf_probe))){
  SNPs<-NULL
  if (dim(info_markers_C2[info_markers_C2$scaf_probe==unique(info_markers_C2$scaf_probe)[i],])[1]==1){
    to_onemap_C2_probe<-rbind(to_onemap_C2_probe,info_markers_C2[info_markers_C2$scaf_probe==unique(info_markers_C2$scaf_probe)[i],])
  }
  else {
    SNPs<- info_markers_C2[info_markers_C2$scaf_probe==unique(info_markers_C2$scaf_probe)[i],]
    SNPs<-SNPs[SNPs$miss==min(SNPs$miss),]
    SNPs<-SNPs[SNPs$maf==max(SNPs$maf),][1,]
    to_onemap_C2_probe<-rbind(to_onemap_C2_probe, SNPs)
  }
}

Onemap_C2_probe<-NULL
Onemap_C2_probe$cross<-NULL
Onemap_C2_probe$marker<-NULL
for (i in 1:dim(to_onemap_C2_probe)[1]){
  Onemap_C2_probe$marker[i]<-paste("*",to_onemap_C2_probe$scaf_probe[i],sep="")
  Onemap_C2_probe$cross[i]<-"D1.11"
}
Onemap_C2_probe$GT<-to_onemap_C2_probe$GT
Onemap_C2_probe<-as.data.frame(Onemap_C2_probe)


##C3
to_onemap_C3_probe<-NULL  
for (i in 1:length(unique(info_markers_C3$scaf_probe))){
  SNPs<-NULL
  if (dim(info_markers_C3[info_markers_C3$scaf_probe==unique(info_markers_C3$scaf_probe)[i],])[1]==1){
    to_onemap_C3_probe<-rbind(to_onemap_C3_probe,info_markers_C3[info_markers_C3$scaf_probe==unique(info_markers_C3$scaf_probe)[i],])
  }
  else {
    SNPs<- info_markers_C3[info_markers_C3$scaf_probe==unique(info_markers_C3$scaf_probe)[i],]
    SNPs<-SNPs[SNPs$miss==min(SNPs$miss),]
    SNPs<-SNPs[SNPs$maf==max(SNPs$maf),][1,]
    to_onemap_C3_probe<-rbind(to_onemap_C3_probe, SNPs)
  }
}

Onemap_C3_probe<-NULL
Onemap_C3_probe$cross<-NULL
Onemap_C3_probe$marker<-NULL
for (i in 1:dim(to_onemap_C3_probe)[1]){
  Onemap_C3_probe$marker[i]<-paste("*",to_onemap_C3_probe$scaf_probe[i],sep="")
  Onemap_C3_probe$cross[i]<-"D1.11"
}
Onemap_C3_probe$GT<-to_onemap_C3_probe$GT
Onemap_C3_probe<-as.data.frame(Onemap_C3_probe)

save.image(".RData")

write.table(Onemap_C1_probe,"Onemap_C1_probe.txt",row.names=F,col.names=F,quote=F)
write.table(Onemap_C2_probe,"Onemap_C2_probe.txt",row.names=F,col.names=F,quote=F)
write.table(Onemap_C3_probe,"Onemap_C3_probe.txt",row.names=F,col.names=F,quote=F)

