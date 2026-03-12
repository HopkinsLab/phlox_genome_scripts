BEGIN { 
    OFS="\t"
    split(rm_samps, rm_array, ",")
    for(i in rm_array){
        check_array[rm_array[i]] = 1
    }

    sp_i = 7
    id_i = 2
}

NR == 1 {
    # Figure out which column to use for pop column based on `pop_col`
    for(i=1; i <= NF; i++){
        if($i == pop_col){
            pop_i = i
            break
        }
    }
}

NR > 1 { if($sp_i == sp && (! check_array[$id_i])){ print $id_i,$pop_i }}
