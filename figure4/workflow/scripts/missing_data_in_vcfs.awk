#!/usr/bin/awk -f
#
# Usage:
#   awk -v wsize=10000 -f missing_windows.awk input.tsv > output.tsv
#
# Input:
#   col1 = chromosome
#   col2 = genomic position (1-based)
#   col3..NF = VCF-style genotypes (./., 0/0, 0/1, 1/1, etc.)
#
# Output:
#   col1 = chromosome
#   col2 = window_start (0-based, inclusive)
#   col3 = window_end   (0-based, exclusive; BED-style)
#   col4.. = number of sites in that window with
#            0, 1, 2, ..., N missing genotypes (N = number of samples)
#

BEGIN {
    OFS = "\t"
}

# First line: determine number of samples (NF - 2)
NR == 1 {
    nSamples = NF - 2
    outstr = "CHROM\tSTART\tEND"
    for(i=0; i <= nSamples; i++){
        outstr = outstr "\tMISSING" i
    }
    printf outstr "\n"
}

{
    chr = $1
    pos = $2 + 0      # 1-based genomic position
    # Convert to 0-based and find window index
    win = int((pos - 1) / wsize)

    # Initialize current window on first record
    if (NR == 1) {
        currentChr = chr
        currentWin = win
    }

    # If we moved to a new window or chromosome, flush previous window
    if (chr != currentChr || win != currentWin) {
        printWindow()
        delete countMissing
        currentChr = chr
        currentWin = win
    }

    # Count number of missing genotypes at this site
    nMiss = 0
    for (i = 3; i <= NF; i++) {
        # Missing if "./." or ".|."
        if ($i ~ /^\.[\/|]\./) {
            nMiss++
        }
    }
    countMissing[nMiss]++
}

END {
    # Flush last window if there was any input
    if (NR > 0) {
        printWindow()
    }
}

# Print current window stats
function printWindow(   start, end, m, val) {
    if (currentChr == "") return  # no data

    start = currentWin * wsize
    end   = start + wsize

    # Chrom, window start, window end
    printf "%s\t%d\t%d", currentChr, start, end

    # Counts for 0..nSamples missing genotypes
    for (m = 0; m <= nSamples; m++) {
        val = (m in countMissing) ? countMissing[m] : 0
        printf "\t%d", val
    }
    printf "\n"
}
