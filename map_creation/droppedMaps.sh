#!/bin/bash -l

#SBATCH -c 4
#SBATCH --mem=50G
#SBATCH -A 
#SBATCH --mail-user 
#SBATCH --mail-type=ALL

set -x

module load R

LG=$1

Rscript DroppedMap.R ${LG}
