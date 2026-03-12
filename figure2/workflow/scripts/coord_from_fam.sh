#!/usr/bin/env bash

metadata=$1
fam=$2
output=$3

latcol=$(  head -n 1 $metadata | tr '\t' '\n' | grep -n '^lat$'  | cut -f 1 -d ':' )
longcol=$( head -n 1 $metadata | tr '\t' '\n' | grep -n '^long$' | cut -f 1 -d ':' )

for samp in `cut -f 2 -d ' ' $fam`; do
    awk -v samp=$samp -v FS="\t" -v OFS=" " -v latcol=$latcol -v longcol=$longcol \
        'NR > 1 && $2 == samp {print $latcol,$longcol}' $metadata >> $output
done
