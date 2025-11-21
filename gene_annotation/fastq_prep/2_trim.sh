#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=4G
#SBATCH -p serial_requeue
#SBATCH -o logs/%x.%A_%a.out
#SBATCH -e logs/%x.%A_%a.err
#SBATCH -t 0-04:00
#SBATCH -J trim
set -e 

i=${SLURM_ARRAY_TASK_ID:-1}
samp=$(sed -n "${i}p" sample_names.list)

fq_1=results/merged_fastqs/${samp}_merged_R1.fastq.gz
fq_2=results/merged_fastqs/${samp}_merged_R2.fastq.gz

conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda activate trim

out_dir=results/trimmed
if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi


trim_galore \
    --paired \
    --retain_unpaired \
    --phred33 \
    --output_dir $out_dir \
    --length 36 \
    -q 5 \
    --stringency 1 \
    -e 0.1 \
    $fq_1 $fq_2
