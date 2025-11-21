#!/usr/bin/env bash

conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda env create -n trim envs/trim.yml

# Merge fastqs
sbatch -a 1-46 1_merge_fastqs.sh 1
sbatch -a 1-46 1_merge_fastqs.sh 2

# Trim
sbatch -a 1-46 2_trim.sh
