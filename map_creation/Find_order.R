load("/mnt/picea/home/cbernhardsson/genetic-maps/Spruce_C1_new/C1_cleaned.RData")

library(BatchMap)

set.map.fun("kosambi")

args <- commandArgs(trailingOnly = TRUE)

message(args[1])

LG_in<-LGs[[ args[1] ]]

LG_rec<-record.parallel(LG_in,times=16,cores=16)

print(LG_rec$seq.num)

size_pred <- pick.batch.sizes(LG_rec, size = 50,overlap = 30, around = 10)

map <- map.overlapping.batches(input.seq=LG_rec,
                              size=size_pred,
                              overlap=30,
                              fun.order=ripple.ord,
                              phase.cores=4,
                              ripple.cores=8,
                              ws=10,
                              max.dist = 25,
                              max.tries = 3,
                              min.tries = 1,
                              method="one",
                              optimize="likelihood",
                              verbosity=c("order","batch"))



save(map,file=paste(args[1],"rec16_10wRippled_Map.RData",sep="_"))