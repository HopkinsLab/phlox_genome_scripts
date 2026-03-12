# Flip alleles from negative to positive strand
function flip_strand(x){
    if(x == "A"){
        x = "T"
    } else if(x == "T"){
        x = "A"
    } else if(x == "C"){
        x = "G"
    } else if(x == "G"){
        x = "C"
    }
    return x
}


BEGIN {
    FS = "\t"
    OFS = "\t"
    CHROM = 1
    POS = 2
    ID = 3
    REF = 4
    ALT = 5
    QUAL = 6
    FILTER = 7
    INFO = 8
    FORMAT = 9
}

/^#/ {
    if($1 == "#CHROM"){
        # Print new INFO line
        print "##INFO=<ID=STRAND_FLIPPED,Number=0,Type=Flag,Description=\"Flipped strand of REF and ALT alleles called by STACKS\">"
    }
    print $0
}

!/^#/ {
    if(sub(/-$/, "+", $ID)){
        if($9 == "GT"){
            $REF = flip_strand($REF)
            $ALT = flip_strand($ALT)
            $INFO = $INFO ";STRAND_FLIPPED"
        }
    }
    print $0
}
