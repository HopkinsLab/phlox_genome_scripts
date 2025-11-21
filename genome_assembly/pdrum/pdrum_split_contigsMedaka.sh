#split file to files with ~4000 contig names in each
awk '{print $1}' curated.with.repeats.5kb.fasta.fai | split -l 4000 - curated.with.repeats.5k

#Replace newline with space
sed -i ':a;N;$!ba;s/\n/ /g' curated.with.repeats.5kba*

#Make file of filenames for array jobs (7 files total)
ls $PWD/*.5kba[a-z] > contigs.fofn

p=$(sed "${SLURM_ARRAY_TASK_ID}q;d" contigs.fofn)
contigs=$(echo $p)