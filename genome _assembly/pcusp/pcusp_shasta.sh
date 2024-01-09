#!/bin/bash
#SBATCH -J shasta
#SBATCH -o out
#SBATCH -e err
#SBATCH -p ultramem
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -t 7-00:00
#SBATCH --mem=1600G

#gunzip -c /n/holyscratch01/informatics/dkhost/phlox_SV/pcusp_Hopkins_801_19_2/pcusp_Hopkins_801_19_2_reads.fastq.gz > pcusp_Hopkins_801_19_2_reads.fastq

/n/holylfs05/LABS/informatics/Users/dkhost/shasta-Linux-0.11.1 --threads 32 --memoryMode anonymous --memoryBacking 4K --config Nanopore-Plants-Apr2021.conf --input pcusp_Hopkins_801_19_2_reads.fastq