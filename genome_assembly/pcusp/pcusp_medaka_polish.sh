#!/bin/bash
#SBATCH -J medaka_array # A single job name for the array
#SBATCH -a 1-26
#SBATCH --partition=gpu # Partition
#SBATCH -N 1 # Number of Nodes required
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G # Memory request
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o medaka_%A_%a.out # Standard output
#SBATCH -e medaka_%A_%a.err # Standard error

module load intel/19.0.5-fasrc01

p=$(sed "${SLURM_ARRAY_TASK_ID}q;d" contigs_fofn.txt)
CONTIGS=$(cat ${p})

BAM="/n/holyscratch01/informatics/dkhost/phlox_cusp_assembly/pcusp_Hopkins_801_reads_vs_pcusp_v0.2.mapped.sorted.bam"

singularity exec --cleanenv --nv /n/holylfs05/LABS/informatics/Users/dkhost/medaka_v1.7.2_gpu.sif medaka consensus $BAM consensus_probs.${SLURM_ARRAY_TASK_ID}.hdf --batch_size 500 --threads 8 --region $CONTIGS --model r941_sup_plant_g610