#!/bin/bash
#SBATCH -n 1
#SBATCH --mem 15000
#SBATCH -p serial_requeue,shared
#SBATCH -J tsIntegrate
#SBATCH -o logs/tdecoder_%x_%A.out
#SBATCH -e logs/tdecoder_%x_%A.err
#SBATCH -t 12:00:00


conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda activate transdecoder_5.5.0

transcriptomefasta=$1
transcriptsgtf=$2 
tdecgff=${transcriptomefasta}.transdecoder.gff3

python IntegrateStringtieNcrna2StringtieTdecGffCdsUtr.py \
    -tdecgff $tdecgff -gtf $transcriptsgtf -gff3out ${tdecgff%.*}.IntegrateNcrna.gff3

