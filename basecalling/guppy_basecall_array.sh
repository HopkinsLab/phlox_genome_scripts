#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH -t 7-00:00                # Runtime in D-HH:MM
#SBATCH --partition=gpu
#SBATCH -a 1-8
#SBATCH --mem=64GB
#SBATCH --gres=gpu:1
#SBATCH -J guppy
#SBATCH -o outfile.%A.out # File to which STDOUT will be written
#SBATCH -e outfile.%A.err # File to which STDERR will be written


module load cuda/10.1.243-fasrc01

GUPPY="/n/holylfs05/LABS/informatics/Users/dkhost/ont-guppy_v6.4.6/bin/guppy_basecaller"

LINE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" SAMPLES_FOFN.txt)
FAST5=$(printf "$LINE" | cut -f2)
OUTDIR=$(printf "$LINE" | cut -f1)

mkdir -p $OUTDIR

ls $FAST5 | $GUPPY --disable_qscore_filtering --config dna_r9.4.1_450bps_sup_plant.cfg --records_per_fastq 0 --num_callers 8 --gpu_runners_per_device 16 -s $OUTDIR --recursive --compress_fastq --device cuda:0