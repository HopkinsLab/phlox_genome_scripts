#!/usr/bin/env bash

conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$conda_setup"
conda env create -n transdecoder_5.5.0 envs/transdecoder.yml
conda env create -n blast_2.12.0 envs/blast.yml

merge_gtf="../hisat2_stringtie/results/stringtie/hisat2/merged_assembly/phlox_stringtie-hisat2_merge.gtf"
transcripts_fa="phlox_stringtie-hisat2_merge.transcripts.fa"

# Step 1
sbatch ./1_TransdecoderBuildTranscriptomeFastaFromGtf.sh $merge_gtf $transcripts_fa

# Step 2
sbatch ./2_Transdecoder_longorfs.sh $transcripts_fa

# Step 3
cd blastp/db
sbatch makedb.sh

cd ..
sbatch ./Tdecblastp2UniprotTrembl_Ericales.sh $transcripts_fa
./collect_parts.sh $transcripts_fa

# Step 4
cd ..
sbatch ./4_Transdecoder_predict.sh $transcripts_fa

# Step 5
sbatch ./5_Transdecoder_CreateGenomeAnnotationFile.sh $transcripts_fa $merge_gtf

sed 's/;Name.*//' ${transcripts_fa}.transdecoder.genome.gff3 > \
 phlox_v1.0.stringtie-hisat2_transdecoder-proteincoding.gff3

# Step 6
sbatch ./6_Transdecoder_IntegrateNcrna.sh $transcripts_fa $merge_gtf
