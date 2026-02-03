#!/usr/bin/env bash
#
#SBATCH -J proc_summary_data
#SBATCH -c 1
#SBATCH --mem-per-cpu=10G
#SBATCH -t 0-00:30
#SBATCH -p serial_requeue
#SBATCH -o logs/%x.%A_%a.o
#SBATCH -e logs/%x.%A_%a.e
set -e

conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda activate popgen

module load gcc/12.2.0-fasrc01
module load R/4.2.2-fasrc01

if [ -d /scratch ]; then
    TMPDIR=/scratch
fi

aggr_method=$( echo ${1:-uniq} | tr '[:upper:]' '[:lower:]' )


if [[ $aggr_method == "pressubsamp" || $aggr_method == "pressub" || $aggr_method == "pressamp" ]]; then
    rep_i=${SLURM_ARRAY_TASK_ID:-1}
    use_seed=${RANDOM}
    
    out_dir="te_sv_counts/presSubsamp"
    if [ ! -d $out_dir ]; then
        mkdir -p $out_dir
    fi

    echo $use_seed > $out_dir/seed.${rep_i}.txt

    ./proc_summary_data.R $aggr_method $out_dir $rep_i $use_seed
else
    out_dir="te_sv_counts"
    if [ ! -d $out_dir ]; then
        mkdir -p $out_dir
    fi

    ./proc_summary_data.R $aggr_method $out_dir
fi

