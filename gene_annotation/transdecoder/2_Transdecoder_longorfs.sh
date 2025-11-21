#!/bin/bash
#SBATCH -n 1
#SBATCH --mem 15000
#SBATCH -p serial_requeue,shared
#SBATCH -o logs/tdecod_longorfs_%A.out
#SBATCH -e logs/tdecod_longorfs_%A.err
#SBATCH -t 23:00:00
#SBATCH -J tdeclongorfs

conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda activate transdecoder_5.5.0


transcriptfasta=$1
echo "params: $transcriptfasta"

TransDecoder.LongOrfs -t $transcriptfasta

# Spilt the protein sequences for downstream blasting
mkdir -p blastp/$transcriptfasta.parts
split -a 3 -l 2000 --numeric-suffixes=1 --additional-suffix=.fasta $transcriptfasta.transdecoder_dir/longest_orfs.pep blastp/$transcriptfasta.parts/longest_orfs.pep.
