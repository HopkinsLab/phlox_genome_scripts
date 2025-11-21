#!/bin/bash
#SBATCH -n 1
#SBATCH --mem 15000
#SBATCH -p serial_requeue,shared
#SBATCH -o logs/tdecoder_ts2genome_%A.out
#SBATCH -e logs/tdecoder_ts2genome_%A.err
#SBATCH -t 16:00:00
#SBATCH -J ts2genome


conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda activate transdecoder_5.5.0

transcriptomefasta=$1
transcriptsgtf=$2 
transcriptsgff3=${transcriptsgtf%.*}.gff3

# Convert transcripts gtf to gff3
gtf_to_alignment_gff3.pl $transcriptsgtf > $transcriptsgff3

cdna_alignment_orf_to_genome_orf.pl \
    ${transcriptomefasta}.transdecoder.gff3 \
    $transcriptsgff3 \
    $transcriptomefasta \
    > ${transcriptomefasta}.transdecoder.genome.gff3

