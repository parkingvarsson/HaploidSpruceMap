read.table("multi_marker_scafs_pre_drop.txt",head=F)->multi_marker
read.table("haploid_scaffolds.txt",head=F)->scaf_length
names(scaf_length)<-c("Scaffold","Length")

read.table("intra_split_scaffolds.txt",head=F)->intra_split
read.table("inter_split_scaffolds.txt",head=F)->inter_split

scaf_length_multi<-scaf_length[scaf_length$Scaffold %in% multi_marker$V1,]
scaf_length_split<-scaf_length[scaf_length$Scaffold %in% inter_split$V1 |scaf_length$Scaffold %in% intra_split$V1,]

scaf_length_multi$type<-"multi"
scaf_length_split$type<-"split"

typeLength<-rbind(scaf_length_multi[,2:3],scaf_length_split[,2:3])

install.packages("ggplot2")
library(ggplot2)

ggplot(typeLength, aes(Length, fill = type)) + geom_density(alpha = 0.2)

plot(density(scaf_length_multi$Length),xlim=c(389,161531),las=1,xlab="",ylab="",main="Scaffold Lengths")
points(density(scaf_length_split$Length),type="l",col="red")
legend("topright", c("Multi marker scaffolds","Split Scaffolds"),col=c("black","red"),pch=19,bty="n",cex=0.7)

