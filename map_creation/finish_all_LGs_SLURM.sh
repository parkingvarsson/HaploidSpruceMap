#!/bin/bash -l

#SBATCH -c 20
#SBATCH --mem=150G
#SBATCH -A 
#SBATCH --mail-user 
#SBATCH --mail-type=ALL

set -x

module load R

LG=$1

BASE_CALL="load('/mnt/picea/projects/spruce/genetic-maps/Spruce_C3_probe_map/.RData');library(onemap);set.map.fun('kosambi');"

CALL="LG_${LG}_rec40<-record.parallel.rfc(LG_${LG},times=40,cores=20);save(LG_${LG}_rec40,file='LG_${LG}_rec40.RData')"

FINAL="Rscript -e \"${BASE_CALL}${CALL}\""

eval $FINAL
