A snakemake (version 7) pipeline to run three long read SV callers (SVIM, Sniffles2 and cuteSV), merge & take intersection of all calls, and force call SVs with Sniffles2.

### Useage
```
snakemake -s Snakefile_v2 --profile profiles/slurm
```
