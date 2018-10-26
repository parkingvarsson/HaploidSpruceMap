#!/bin/bash -l

#SBATCH -p node
#SBATCH -n 16
#SBATCH -A u2017005
#SBATCH --mail-user carolina.bernhardsson@umu.se
#SBATCH --mail-type=ALL
#SBATCH --mem 128G

set -x

module load R

LG=$1

Rscript Find_order.R ${LG}
