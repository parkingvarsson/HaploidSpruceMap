library(onemap)

read_outcross_cpp("Onemap_C3_probe.txt")->C3.out
save.image(".RData")
twopts_C3 <- rf.2pts(C3.out, LOD=8, max.rf=0.35)

mark.all_C3 <- make.seq(twopts_C3, "all")
save.image(".RData")

class(C3.out)[2]<- "outcross"
class(twopts_C3)[2] <- "outcross"

LGs_C3<- group(mark.all_C3)
print(LGs_C3,detailed=F) ## MAKE SURE THAT THE DATA GROUPS INTO 12 LGs, OTHERWISE CHANGE LOD SCORE within group() command

## make LGs into sequences
for (i in 1:LGs_C3$n.groups){
  print(i)
  assign(paste("LG",i,sep="_"),make.seq(LGs_C3,i))
}

set.map.fun("kosambi")

save.image(".RData")
