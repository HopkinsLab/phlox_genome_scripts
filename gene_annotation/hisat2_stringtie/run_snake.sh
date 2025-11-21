#!/usr/bin/env bash
#
#SBATCH -J runsnake
#SBATCH -c 1
#SBATCH --mem-per-cpu=4G
#SBATCH -t 4320
#SBATCH -p shared

source ~/.bashrc 
conda activate snakemake

snakemake --profile slurm --snakefile Snakefile_RNAseq
