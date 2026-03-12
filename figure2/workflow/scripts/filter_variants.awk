BEGIN {
    FS = "\t"
    OFS = "\t"

    # Get sample populations
    while((getline line < popmap) > 0 ){
        split(line, arr, "\t")    
        pop_info[arr[1]] = arr[2]
        nsamps[arr[2]] += 1
    } 
    close(popmap)
}

/^#/ { 
    print $0
    if($1 == "#CHROM"){ # Get sample names
        for(s=10; s <= NF; s++){
            sampnames[s] = $s
        }
        nsamps["all"] = NF - 9
    }
}

!/^#/ {
    keep_record = 1
    

    # Initialize missing counters
    for(pop in nsamps){
        nmiss[pop] = 0
    }

    if($9 == "GT"){
        # for(s=10; s <= NF; s++){
        #     is_miss = ($s == "./.")
        #     samp = sampnames[s]
        #     pop  = pop_info[samp]
        #     if(is_miss){ # Count missing genotype
        #         nmiss["all"] += 1
        #         nmiss[pop] += 1
        #     }
        # }
    } else {
        split($9, fmt_names, ":")
        for(s=10; s <= NF; s++){
            is_miss = (substr($s, 1, 3) == "./.")
            samp = sampnames[s]
            pop  = pop_info[samp]
            if(is_miss){ # Count missing genotype
                nmiss["all"] += 1
                nmiss[pop] += 1
            } else if($s != "0/0"){ # Filter out bad genotypes
                split($s, fmt_stats, ":")
                for(i in fmt_names){
                    if(fmt_names[i] == "DP" && (fmt_stats[i] < min_DP || fmt_stats[i] > max_DP )){
                        $s = "./."
                        nmiss["all"] += 1
                        nmiss[pop] +=1
                        break
                    } else if(fmt_names[i] == "GQ" && (fmt_stats[i] < min_GQ)){
                        $s = "./."
                        nmiss["all"] += 1
                        nmiss[pop] +=1
                        break
                    }
                }
            }
        }
        # Check if we should filter out record entirely
        for(pop in nsamps){
            if(pop == "all"){
                cutoff = min_samp_overall
            } else {
                cutoff = min_samp_per_pop
            }
            # print pop":"nmiss[pop]","nsamps[pop]
            if(nmiss[pop] / nsamps[pop] > 1 - cutoff){
                keep_record = 0
                break
            }
        }
    }


    if(keep_record){
        print $0
    }
}
