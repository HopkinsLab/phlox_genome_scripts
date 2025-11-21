#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH --mem 1G
#SBATCH -p serial_requeue
#SBATCH -t 00:00:30
#SBATCH -J collect

transcriptfa=$1
out_dir="${transcriptfa}.parts"
out_fn="${out_dir}/longorfs.blastp2SwissProt.concat.tsv"

if [[ ! -s ${out_fn} ]]; then
    for i in {1..273}; do 
        cat ${out_dir}/longorfs.blastp2SwissProt.$(printf %03d $i).tsv >> $out_fn
    done
fi
