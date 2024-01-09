#!/bin/bash
#SBATCH -n 10                # Number of cores
#SBATCH -t 14-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p holy-cow   # Partition to submit to
#SBATCH --mem=120G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o outfile.out # File to which STDOUT will be written
#SBATCH -e outfile.err # File to which STDERR will be written

export SLURM_TMPDIR=/n/holyscratch01/informatics/dkhost/tmp
export TMPDIR=/n/holyscratch01/informatics/dkhost/tmp

source ~/mambaforge/bin/activate minimap
REF="pcusp_ragtag_vs_pdrum/ragtag.scaffold.fasta"
SAMPLE="pcusp_Hopkins_801"

minimap2 -ax map-ont --secondary=no --sam-hit-only -K 6G -I 32G -t 10 $REF pcusp_Hopkins_801_19_2_reads.fastq > ${SAMPLE}_reads_vs_pcusp_v0.2.sam

source ~/mambaforge/bin/activate samtools

samtools view -T $REF -h -b -F 4 ${SAMPLE}_reads_vs_pcusp_v0.2.sam | samtools sort -@ 10 -T /n/holyscratch01/informatics/dkhost/tmp -o ${SAMPLE}_reads_vs_pcusp_v0.2.mapped.sorted.bam -
samtools index -c ${SAMPLE}_reads_vs_pcusp_v0.2.mapped.sorted.bam