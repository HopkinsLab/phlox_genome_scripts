#!/usr/bin/env bash

i=$( echo "$1" | md5sum | cut -f 1 -d ' ' | sed 's/[a-z]//g' )

# Truncate to last 5 digits
if [[ ${#i} -gt 5 ]]; then
    i=${i: -5}
fi

echo $i
