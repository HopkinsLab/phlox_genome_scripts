#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=4G
#SBATCH -p serial_requeue
#SBATCH -o logs/%x.%A_%a.out
#SBATCH -e logs/%x.%A_%a.err
#SBATCH -t 0-01:00
#SBATCH -J merge_fastqs
set -e 

i=${SLURM_ARRAY_TASK_ID:-1}
samp=$(sed -n "${i}p" sample_names.list)

pair=R${1:-1}

fq_arr=( $(awk -v samp=$samp -v pair=$pair '$1==samp $2==pair {print $3}' raw_fastqs.paths.tsv) )

out_dir=results/merged_fastqs
if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

if [ ${#fq_arr[@]} -eq 1 ]; then
    ln -s ${fq_arr[0]} ${out_dir}/${samp}_merged_${pair}.fastq.gz
else
    cat ${fq_arr[@]} > ${out_dir}/${samp}_merged_${pair}.fastq.gz
fi
