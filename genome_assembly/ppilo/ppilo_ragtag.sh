#!/bin/bash
#SBATCH -J ragtag_array # A single job name for the array
#SBATCH --partition=holy-cow # Partition
#SBATCH -n 12 # Number of Nodes required
#SBATCH --mem=200G # Memory request
#SBATCH -t 30-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o ragtag.out # Standard output
#SBATCH -e ragtag.err # Standard error

source ~/mambaforge/bin/activate ragtag

ASSEM="ShastaRun/Assembly.fasta"
REF="/n/holylfs05/LABS/hopkins_lab/Lab/Phlox_resources/WG_assemblies/assemblies/drummondii/phlox_flye.v1.0.FINAL.fasta"

ragtag.py scaffold --remove-small -f 10000 -w -o ppilo_ragtag_vs_pdrum -t 12 -u --mm2-params='-x asm10 -I 20G' $REF $ASSEM