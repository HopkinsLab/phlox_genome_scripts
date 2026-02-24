# Repeat figures

Code for reproducing repeat figures (largely Figure 3 in initial draft manuscript)

## Dependencies

### R

R and Rmarkdown packages were executed using R version 4.2.2.

Required packages:
* `data.table`
* `string`
* `ggpubr`
* `patchwork`

## Execution

First, compile te-sv counts:

```
sbatch proc_summary_data.sh 1sp 
sbatch proc_summary_data.sh 2sp
sbatch proc_summary_data.sh 3sp
sbatch proc_summary_data.sh pres
sbatch -a 1-500 proc_summary_data.sh presSubsamp
```

Next, plotting. This is all done in the Rmarkdown file so refer there.

```
../utils/render_rmd.R figure3.Rmd
```
