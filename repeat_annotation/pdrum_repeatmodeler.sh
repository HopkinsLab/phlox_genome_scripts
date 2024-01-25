#!/bin/bash
#SBATCH -n 20                # Number of cores
#SBATCH -N 1
#SBATCH -t 14-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p holy-cow   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o repmod.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e repmod.err  # File to which STDERR will be written, %j inserts jobid


/n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools.sh --container /n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools-latest.sif -- BuildDatabase -name phlox_flye.v1.0.FINAL.lg_split /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_resources/WG_assemblies/assemblies/drummondii/phlox_flye.v1.0.FINAL.lg_split.fasta

/n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools.sh --container /n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools-latest.sif -- RepeatModeler -database phlox_flye.v1.0.FINAL.lg_split -pa 20 -LTRStruct