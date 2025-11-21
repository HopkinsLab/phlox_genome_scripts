#!/bin/bash
#SBATCH -J blastp
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p serial_requeue,shared # Partition 
#SBATCH -t 02:00:00
#SBATCH --array=1-216
#SBATCH --mem=8000
#SBATCH -o logs/blstp_%A_%a.out 
#SBATCH -e logs/blstp_%A_%a.err 

transcriptfasta=$1

db_faa="db/Ericales_uniprotTrembl_2022.10.10.faa"
out_dir="${transcriptfasta}.parts"

conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda activate blast_2.12.0 

i=$( printf %03d ${SLURM_ARRAY_TASK_ID} )

blastp -max_target_seqs 5 -num_threads 16  -evalue 1e-4 \
    -query ${out_dir}/longest_orfs.pep.${i}.fasta -outfmt 6 \
    -db $db_faa > ${out_dir}/longorfs.blastp2SwissProt.${i}.tsv

