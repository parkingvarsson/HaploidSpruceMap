read.table("../Consensus_maps.txt",head=T)->maps
read.table("haploid_scaffolds.txt",head=F)->scaf_length
names(scaf_length)<-c("Scaffold","Length")

read.table("intra_split_scaffolds.txt",head=F)->intra_split
read.table("inter_split_scaffolds.txt",head=F)->inter_split

###Scaffolds with splits of 5-10 cM in consensus map or where markers on each side of split are from different framework maps 
maps[maps$Scaffold %in% c("MA_9458","MA_281725","MA_10431182","MA_10431315","MA_10432328","MA_13639"),]
-----
#  LG        Marker Consensus      C3      C2      C1    Scaffold Probe
#1      1 MA_10431182:1     0.000   0.000      NA      NA MA_10431182     1   NOT A TRUE SPLIT
#9      1 MA_10431182:2    29.917      NA   2.379      NA MA_10431182     2
#2598   2 MA_10431315:2    63.012  69.218      NA  43.591 MA_10431315     2   NOT A TRUE SPLIT
#2599   2 MA_10431315:4    63.012  69.218      NA      NA MA_10431315     4
#2605   2 MA_10431315:1    63.442  69.218      NA  44.058 MA_10431315     1
#2607   2 MA_10431315:3    70.104  69.448      NA      NA MA_10431315     3
#4965   3 MA_10432328:4   180.242 208.290 212.489 175.669 MA_10432328     4   NOT A TRUE SPLIT
#4975   3 MA_10432328:2   185.646 207.721      NA 174.135 MA_10432328     2
#4978   3 MA_10432328:5   185.892 206.938      NA      NA MA_10432328     5
#4980   3 MA_10432328:3   186.931 208.290      NA 176.146 MA_10432328     3
#8004   5    MA_13639:1    98.503 111.164 116.049  96.516    MA_13639     1   TRUE SPLIT!
#8009   5    MA_13639:2   104.242 116.904 124.725      NA    MA_13639     2
#15420  9     MA_9458:2   195.663 188.926      NA      NA     MA_9458     2   NOT A TRUE SPLIT
#15425  9     MA_9458:1   201.239 188.926 201.459      NA     MA_9458     1
#16028 10   MA_281725:2     0.000   0.000      NA      NA   MA_281725     2   NOT A TRUE SPLIT
#16061 10   MA_281725:1    39.678      NA      NA   3.902   MA_281725     1
----

#install.packages("gdata")
library(gdata)
intra_split<-intra_split[!(intra_split$V1 %in% c("MA_9458","MA_281725","MA_10431182","MA_10431315","MA_10432328")),]
intra_split<-drop.levels(intra_split)  
  
multi_marker_scafs<-NULL
for (i in 1:length(levels(maps$Scaffold))){
  tmp <- maps[maps$Scaffold == levels(maps$Scaffold)[i],]
  if (dim(tmp)[1] > 1){
    multi_marker_scafs<-rbind(multi_marker_scafs,tmp)
  }
}

multi_marker_scafs<-drop.levels(multi_marker_scafs)



scaf_length_multi<-scaf_length[as.character(scaf_length$Scaffold) %in% as.character(multi_marker_scafs$Scaffold),]
scaf_length_split<-scaf_length[as.character(scaf_length$Scaffold) %in% as.character(inter_split$V1) |as.character(scaf_length$Scaffold) %in% as.character(intra_split),]

scaf_length_multi$type<-"multi"
scaf_length_split$type<-"split"

pdf("Length_multiMarker_vs_Split.pdf")
boxplot(scaf_length_multi$Length/1e3,scaf_length_split$Length/1e3,notch = T,las=1,ylab="Scaffold length (kb)",col=c(rgb(0.3,0.3,0.3),rgb(0.7,0.7,0.7)))
df<-t.test(scaf_length_multi$Length,scaf_length_split$Length)
#legend("topright",c("N = 4859", "N = 185"),pch =19, col=c(rgb(0.3,0.3,0.3),rgb(0.7,0.7,0.7)),bty="n")
mtext(side=1,at=c(1,2),c("All multi-marker scaffolds","Split scaffolds"),line=1)
dev.off()
rm(df)

save.image("MultiScaf_vs_SplitScaf.RData")