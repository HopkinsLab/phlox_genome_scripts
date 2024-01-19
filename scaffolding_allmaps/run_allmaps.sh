#!/bin/bash

#SBATCH -t 7-00:00:00
#SBATCH --mem=64G
#SBATCH -n 8
#SBATCH --partition holy-cow
#SBATCH -J allmaps

source /n/holyscratch01/informatics/dkhost/holylfs_shortcut/phlox_scaffolding/allmaps/venv/bin/activate

python -m jcvi.assembly.allmaps mergebed GeneticMap.bed -o GeneticMap_merged.bed
python -m jcvi.assembly.allmaps path GeneticMap_merged.bed /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_assembly_collab/pdrum_ref/phlox_flye_v5.1.bionano.all_tigs.scaff33623_split.fasta