##Import data
read.table("Cluster1-all_rec40_maps.txt",head=F)->C1_all_rec40_maps
read.table("Cluster2-all_rec40_maps.txt",head=F)->C2_all_rec40_maps
read.table("Cluster3_40rec-Maps.txt",head=F)->C3_all_rec40_maps
names(C1_all_rec40_maps)<-c("LG","Marker","Size")
names(C2_all_rec40_maps)<-c("LG","Marker","Size")
names(C3_all_rec40_maps)<-c("LG","Marker","Size")
C1_all_rec40_maps$scaffold<- as.factor(t(as.data.frame(strsplit(as.character(C1_all_rec40_maps$Marker),split=":")))[,1])
C1_all_rec40_maps$probe<- t(as.data.frame(strsplit(as.character(C1_all_rec40_maps$Marker),split=":")))[,2]
C2_all_rec40_maps$scaffold<- as.factor(t(as.data.frame(strsplit(as.character(C2_all_rec40_maps$Marker),split=":")))[,1])
C2_all_rec40_maps$probe<- t(as.data.frame(strsplit(as.character(C2_all_rec40_maps$Marker),split=":")))[,2]
C3_all_rec40_maps$scaffold<- as.factor(t(as.data.frame(strsplit(as.character(C3_all_rec40_maps$Marker),split=":")))[,1])
C3_all_rec40_maps$probe<- t(as.data.frame(strsplit(as.character(C3_all_rec40_maps$Marker),split=":")))[,2]

##split maps into LGs
C1_LG1<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 1,]
C1_LG2<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 2,]
C1_LG3<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 3,]
C1_LG4<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 4,]
C1_LG5<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 5,]
C1_LG6<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 6,]
C1_LG7<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 7,]
C1_LG8<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 8,]
C1_LG9<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 9,]
C1_LG10<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 10,]
C1_LG11<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 11,]
C1_LG12<-C1_all_rec40_maps[C1_all_rec40_maps$LG == 12,]

library(gdata)
C1_LG1<-drop.levels(C1_LG1)

##non-unique scaffolds per LG
splitted_C1_scaffolds<-NULL
for(i in 1:11){
    for(j in 2:12){
      splitted_C1_scaffolds<-c(splitted_C1_scaffolds,as.character(unique(paste("C1_LG",i,"$scaffold[C1_LG",i,"$scaffold %in% C1_LG",j,"$scaffold]",sep=""))))
    }
}
splitted_C1_scaffolds<-as.factor(splitted_C1_scaffolds)

#Which LGs are correlated between clusters
for (i in 1:12){
  for(j in 1:12){
    print(c("C1",i,"C2",j,length(as.character(unique(paste("C1_LG",i,"$scaffold[C1_LG",i,"$scaffold %in% C2_LG",j,"$scaffold]",sep=""))))))
    print(c("C1",i,"C3",j,length(as.character(unique(paste("C1_LG",i,"$scaffold[C1_LG",i,"$scaffold %in% C3_LG",j,"$scaffold]",sep=""))))))
    print(c("C2",i,"C3",j,length(as.character(unique(paste("C2_LG",i,"$scaffold[C2_LG",i,"$scaffold %in% C3_LG",j,"$scaffold]",sep=""))))))
  }
}

##List all markers from scaffolds splitted over LGs for all clusters
split_scaf_markers_C3<-NULL
for (i in 1:length(levels(C3_all_rec40_maps$scaffold))){
  if (levels(C3_all_rec40_maps$scaffold)[i] %in% splitted_C3_scaffolds){
    split_scaf_markers_C3<-rbind(split_scaf_markers_C3,C3_all_rec40_maps[C3_all_rec40_maps$scaffold == levels(C3_all_rec40_maps$scaffold)[i],])
  }
}
split_scaf_markers_C3<-as.data.frame(split_scaf_markers_C3)
split_scaf_markers_C3<-drop.levels(split_scaf_markers_C3)

#Scaffolds with > 1 marker per LG
C3_common_markers<-NULL
for (i in 1:length(levels(C3_all_rec40_maps$scaffold))){
  for (j in 1:12){
    if (dim(C3_all_rec40_maps[C3_all_rec40_maps$LG ==j & C3_all_rec40_maps$scaffold==levels(C3_all_rec40_maps$scaffold)[i],])[1] >1){
      C3_common_markers<-rbind(C3_common_markers,C3_all_rec40_maps[C3_all_rec40_maps$LG ==j & C3_all_rec40_maps$scaffold==levels(C3_all_rec40_maps$scaffold)[i],])
    }
  }  
}  
C3_common_markers<-drop.levels(C3_common_markers)  
  

#scaffolds with markers >10 cM apart on an LG
suspicious_scaf_C3<-NULL
for (i in 1:length(levels(C3_common_markers$scaffold))){
  for (j in 1:12){
    if (length(C3_common_markers$Size[C3_common_markers$LG == j & C3_common_markers$scaffold == levels(C3_common_markers$scaffold)[i]])>1){
      df<-C3_common_markers[C3_common_markers$LG == j & C3_common_markers$scaffold == levels(C3_common_markers$scaffold)[i],]
      if (max(df$Size) -min(df$Size) > 15){
        suspicious_scaf_C3<-rbind(suspicious_scaf_C3,df)
      }
    }
  }  
}
suspicious_scaf_C3<-drop.levels(suspicious_scaf_C3)
suspicious_scaf_C3$Cluster<-"C3"

all_suspicious_scaffolds<-rbind(suspicious_scaf_C1,suspicious_scaf_C2,suspicious_scaf_C3)
all_suspicious_scaffolds<-drop.levels(all_suspicious_scaffolds)

#which markers are in which cluster
all_markers<-NULL
all_markers$Marker<-unique(c(levels(C1_all_rec40_maps$Marker),levels(C2_all_rec40_maps$Marker),levels(C3_all_rec40_maps$Marker)))
all_markers$C1<-NULL
all_markers$C2<-NULL
all_markers$C3<-NULL
all_markers$Clusters<-0
all_markers<-as.data.frame(all_markers)
for (i in 1:dim(all_markers)[1]){
  if (all_markers$Marker[i] %in% levels(C1_all_rec40_maps$Marker)){
    all_markers$C1[i]<-"YES"
    all_markers$Clusters[i]<-all_markers$Clusters[i]+1
  } else {all_markers$C1[i]<-"NO"}
  if (all_markers$Marker[i] %in% levels(C2_all_rec40_maps$Marker)){
    all_markers$C2[i]<-"YES"
    all_markers$Clusters[i]<-all_markers$Clusters[i]+1
  }else {all_markers$C2[i]<-"NO"}
  if (all_markers$Marker[i] %in% levels(C3_all_rec40_maps$Marker)){
    all_markers$C3[i]<-"YES"
    all_markers$Clusters[i]<-all_markers$Clusters[i]+1
  }else {all_markers$C3[i]<-"NO"}
}

all_markers$scaffold<- as.factor(t(as.data.frame(strsplit(as.character(all_markers$Marker),split=":")))[,1])
all_markers$probe<- t(as.data.frame(strsplit(as.character(all_markers$Marker),split=":")))[,2]



##Drop_markers from LGs
C3_LG12_drop<-NULL
for (i in 1:length(levels(C3_LG12$scaffold))){
  df<-C3_LG12[C3_LG12$scaffold == levels(C3_LG12$scaffold)[i],]
  if (dim(df)[1]>1 & (df$Size[dim(df)[1]]-df$Size[1])<15){
    tmp<-drop.levels(all_markers[all_markers$scaffold == levels(C3_LG12$scaffold)[i] & all_markers$probe %in% df$probe,])
    keep<-as.character(tmp$Marker[tmp$Clusters==max(tmp$Clusters)][1])
    C3_LG12_drop<-c(C3_LG12_drop,as.character(tmp$Marker[tmp$Marker!=keep]))
  }
}


##compare orders between clusters

C1LG1_C2LG11<-NULL
C1LG1_C2LG11$C1<-which(C1_LG1$Marker %in% C2_LG11$Marker)
C1LG1_C2LG11$C2<-NULL
for(i in 1:length(C1LG1_C2LG11$C1)){
  C1LG1_C2LG11$C2[i]<-which(as.character(C2_LG11$Marker) == as.character(C1_LG1$Marker[C1LG1_C2LG11$C1[i]]))
}
C1LG1_C3LG10<-NULL
C1LG1_C3LG10$C1<-which(C1_LG1$Marker %in% C3_LG10$Marker)
C1LG1_C3LG10$C3<-NULL
for(i in 1:length(C1LG1_C3LG10$C1)){
  C1LG1_C3LG10$C3[i]<-which(as.character(C3_LG10$Marker) == as.character(C1_LG1$Marker[C1LG1_C3LG10$C1[i]]))
}
C2LG11_C3LG10<-NULL
C2LG11_C3LG10$C2<-which(C2_LG11$Marker %in% C3_LG10$Marker)
C2LG11_C3LG10$C3<-NULL
for(i in 1:length(C2LG11_C3LG10$C2)){
  C2LG11_C3LG10$C3[i]<-which(as.character(C3_LG10$Marker) == as.character(C2_LG11$Marker[C2LG11_C3LG10$C2[i]]))
}

## Plot heatmap of cluster comparison

pheatmap(comp_matrix,cluster_rows=FALSE,cluster_cols=FALSE,main = "Shared scaffolds between Cluster LGs",gaps_row = c(12,24),gaps_col = c(12,24),display_numbers = matrix(ifelse(m<1,"",m),nrow(m)))
