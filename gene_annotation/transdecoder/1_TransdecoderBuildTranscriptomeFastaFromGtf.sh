#!/bin/bash
#SBATCH -n 1
#SBATCH --mem 25000
#SBATCH -p serial_requeue,shared
#SBATCH -o logs/tdec_gtf2transcriptfasta_%A.out
#SBATCH -e logs/tdec_gtf2transcriptfasta_%A.err
#SBATCH -t 06:00:00
#SBATCH -J gtf2fasta

conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda activate transdecoder_5.5.0

gtfin=$1
tsfastaout=$2
genomefasta="../genomes/phlox_flye.v1.0.FINAL.fasta"

echo "args: $gtfin $genomefasta $tsfastaout"
gtf_genome_to_cdna_fasta.pl $gtfin $genomefasta > $tsfastaout 
