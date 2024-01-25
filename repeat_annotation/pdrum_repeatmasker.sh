#!/bin/bash
#SBATCH -n 20                # Number of cores
#SBATCH -N 1
#SBATCH -t 14-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p holy-cow   # Partition to submit to
#SBATCH --mem=72G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o repmask.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e repmask.err  # File to which STDERR will be written, %j inserts jobid


/n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools.sh --container /n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools-latest.sif -- RepeatMasker -xsmall -pa 9 -gff -no_is -lib ../phlox_flye.v1.0_combined_repeats.fa ../phlox_flye.v1.0.FINAL.fasta