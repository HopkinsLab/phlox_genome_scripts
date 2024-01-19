
#Make file of snps
cut -f 1,2 /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_assembly_collab/linkage_map/res_lepmap/populations.snps.filter_11.call | awk '(NR>=7)' > snps.txt

#Format each linkage group 1..7
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_assembly_collab/linkage_map/res_lepmap/populations.snps.filter_11.sepchrom_lodLim50_t0.0013.js.order.lg1.txt > order1_mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_assembly_collab/linkage_map/res_lepmap/populations.snps.filter_11.sepchrom_lodLim50_t0.0013.js.order.lg2.txt > order2_mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_assembly_collab/linkage_map/res_lepmap/populations.snps.filter_11.sepchrom_lodLim50_t0.0013.js.order.lg3.txt > order3_mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_assembly_collab/linkage_map/res_lepmap/populations.snps.filter_11.sepchrom_lodLim50_t0.0013.js.order.lg4.txt > order4_mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_assembly_collab/linkage_map/res_lepmap/populations.snps.filter_11.sepchrom_lodLim50_t0.0013.js.order.lg5.txt > order5_mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_assembly_collab/linkage_map/res_lepmap/populations.snps.filter_11.sepchrom_lodLim50_t0.0013.js.order.lg6.txt > order6_mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt /n/holylfs05/LABS/hopkins_lab/Lab/Phlox_assembly_collab/linkage_map/res_lepmap/populations.snps.filter_11.sepchrom_lodLim50_t0.0013.js.order.lg7.txt > order7_mapped.txt

#Reformat and combine to single map file
awk 'BEGIN{FS="\t"; OFS="\t"} !/^#/ {print $1,$2,$2+1, "lg1:"$3}' order1_mapped.txt >> GeneticMap.bed
awk 'BEGIN{FS="\t"; OFS="\t"} !/^#/ {print $1,$2,$2+1, "lg2:"$3}' order2_mapped.txt >> GeneticMap.bed
awk 'BEGIN{FS="\t"; OFS="\t"} !/^#/ {print $1,$2,$2+1, "lg3:"$3}' order3_mapped.txt >> GeneticMap.bed
awk 'BEGIN{FS="\t"; OFS="\t"} !/^#/ {print $1,$2,$2+1, "lg4:"$3}' order4_mapped.txt >> GeneticMap.bed
awk 'BEGIN{FS="\t"; OFS="\t"} !/^#/ {print $1,$2,$2+1, "lg5:"$3}' order5_mapped.txt >> GeneticMap.bed
awk 'BEGIN{FS="\t"; OFS="\t"} !/^#/ {print $1,$2,$2+1, "lg6:"$3}' order6_mapped.txt >> GeneticMap.bed
awk 'BEGIN{FS="\t"; OFS="\t"} !/^#/ {print $1,$2,$2+1, "lg7:"$3}' order7_mapped.txt >> GeneticMap.bed