#!/bin/bash
#SBATCH -A snic2017-1-499
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 4-00:00:00
#SBATCH --mail-user carolina.bernhardsson@umu.se
#SBATCH --mail-type=ALL


Chr=$1
hapflk --bfile SwAsp_94samples.south_north_only --chr $Chr -K 15 --nfit=15 --ncpu=5 -p ${Chr}_K15_fit15
