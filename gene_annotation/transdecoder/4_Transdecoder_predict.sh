#!/bin/bash
#SBATCH -c 1
#SBATCH --mem 15000
#SBATCH -p serial_requeue,shared
#SBATCH -o logs/tdecoder_predict_%A.out
#SBATCH -e logs/tdecoder_predict_%A.err
#SBATCH -t 16:00:00
#SBATCH -J tdpred


conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda activate transdecoder_5.5.0

transcriptfasta=$1
blastphits=blastp/${transcriptfasta}.parts/longorfs.blastp2SwissProt.concat.tsv

echo "params: $transcriptfasta"
echo "        $blastphits"

TransDecoder.Predict -t $transcriptfasta --single_best_only --retain_blastp_hits $blastphits

# Outputs files: phlox_stringtie-hisat2_merge.cdna.fa.transdecoder.*
