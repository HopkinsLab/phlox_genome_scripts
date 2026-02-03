# Repeat figures

Code for reproducing repeat figures (largely Figure 3 in initial draft manuscript)

## Dependencies

### R

R and Rmarkdown packages were executed using R version 4.2.2.

Required packages:
* `data.table`
* `string`

## Execution

First, compile te-sv counts:

```
sbatch proc_summary_data.sh 1sp 
sbatch proc_summary_data.sh 2sp
sbatch proc_summary_data.sh 3sp
sbatch proc_summary_data.sh pres
sbatch proc_summary_data.sh presSubsamp
```

Next, plotting.

