bam="/n/holyscratch01/informatics/dkhost/phlox_flye/phlox_flye.medaka/calls_to_draft.bam"
genome="/scratch/phlox_nanopore/phlox_flye/assembly.fasta"

#purge_haplotigs hist -b $bam -g $genome -t 8

# Manually inspected histogram

purge_haplotigs cov -i calls_to_draft.bam.gencov -l 7 -m 32 -h 70

purge_haplotigs purge -g $genome -b $bam -c coverage_stats.csv -t 8 -d