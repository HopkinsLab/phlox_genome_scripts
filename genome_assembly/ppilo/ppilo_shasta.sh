#!/bin/bash
#SBATCH -J shasta
#SBATCH -o out
#SBATCH -e err
#SBATCH -p ultramem,holy-cow
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -t 7-00:00
#SBATCH --mem=1500G

READS1="/n/ngsdata/Nanopore/230302_Hopkins/1704-9/*/fastq*/*.fastq.gz"
READS2="/n/holylfs05/LABS/informatics/Users/dkhost/PHLOX/nanopore_guppy/Hopkins_1704/*.fastq.gz"

zcat $READS1 $READS2 > ppilo_Hopkins_1704_reads.fastq

/n/holylfs05/LABS/informatics/Users/dkhost/shasta-Linux-0.11.1 --threads 32 --memoryMode anonymous --memoryBacking 4K --config Nanopore-Plants-Apr2021.conf --input ppilo_Hopkins_1704_reads.fastq