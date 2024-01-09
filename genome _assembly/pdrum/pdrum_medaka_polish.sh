#!/bin/bash
#SBATCH -J medaka_array # A single job name for the array
#SBATCH -a 1-7
#SBATCH --partition=gpu # Partition
#SBATCH -N 1 # Number of Nodes required
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G # Memory request
#SBATCH -t 5-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o medaka_%A_%a.out # Standard output
#SBATCH -e medaka_%A_%a.err # Standard error

module load intel/19.0.5-fasrc01
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

reads="/n/boslfs02/LABS/informatics/dkhost/phlox_wholegenome/phlox_wg_COMBO.112320.fastq.gz"
assem="/n/holyscratch01/informatics/dkhost/phlox_flye_purged/curated_w_repeats/curated.with.repeats.fasta"

p=$(sed "${SLURM_ARRAY_TASK_ID}q;d" contigs.fofn)
contigs=$(cat ${p})

singularity exec --cleanenv --nv /n/singularity_images/informatics/medaka/medaka_1.2.6-gpu.sif medaka consensus ../phlox_flye.cwr/calls_to_draft.bam ../phlox_flye.cwr/consensus_probs.hdf --model r941_prom_high_g360 --batch_size 100 --threads 8 --region $contigs