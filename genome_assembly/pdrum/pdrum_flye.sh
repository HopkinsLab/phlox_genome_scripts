reads1="/n/boslfs02/LABS/informatics/dkhost/phlox_wholegenome/phlox_wg_COMBO.112320.fastq.gz"

#--nano-raw: indicates uncorrected Nanopore reads. Options for corrected reads.
#-o: path where assembly will be output
#--threads: number of threads used for assembly
time flye --nano-raw $reads1 -g 6.6g -o /scratch/phlox_nanopore/phlox_flye --threads 30