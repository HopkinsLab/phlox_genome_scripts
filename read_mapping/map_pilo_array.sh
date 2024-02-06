#!/bin/bash
#SBATCH -J pilomap # A single job name for the array
#SBATCH -a 1-12
#SBATCH --partition=holy-cow,bigmem_intermediate # Partition
#SBATCH -n 16 # Number of Nodes required
#SBATCH -N 1
#SBATCH --mem=200G # Memory request 
#SBATCH -t 14-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o pilomap_%A_%a.out # Standard output
#SBATCH -e pilomap_%A_%a.err # Standard error


source ~/mambaforge/bin/activate minimap

READDIR=$(sed "${SLURM_ARRAY_TASK_ID}q;d" samples_fastq_fofn.txt)
SAMPLE=$(echo $READDIR | awk 'BEGIN{FS="/"} {print $9}')

export SLURM_TMPDIR=/n/holyscratch01/informatics/dkhost/tmp
export TMPDIR=/n/holyscratch01/informatics/dkhost/tmp

mkdir -p $SAMPLE
cd $SAMPLE

ref="/n/holylfs05/LABS/hopkins_lab/Lab/Phlox_resources/WG_assemblies/assemblies/pilosa/phlox_pilo.v1.fasta"

zcat $READDIR/*.fastq.gz | gzip -1 > ${SAMPLE}_reads_combo.fastq.gz 

minimap2 -ax map-ont --secondary=no -I 30g -K 1G -L -t 15 $ref ${SAMPLE}_reads_combo.fastq.gz > ${SAMPLE}_reads_combo_vs_ppilo_v1.sam