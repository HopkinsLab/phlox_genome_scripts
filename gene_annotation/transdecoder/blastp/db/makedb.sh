#!/bin/bash
#SBATCH -n 1
#SBATCH --mem 15000
#SBATCH -p serial_requeue
#SBATCH -J makedb
#SBATCH -o logs/%x_%A.out
#SBATCH -e logs/%x_%A.err
#SBATCH -t 0-04:00

conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda activate blast_2.12.0 

gunzip Ericales_uniprotTrembl_2022.10.10.faa.gz

makeblastdb -in Ericales_uniprotTrembl_2022.10.10.faa -input_type fasta -dbtype prot 
