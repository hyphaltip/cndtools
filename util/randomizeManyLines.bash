#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 input > output" 1>&2
    exit 1
fi

# Make a temporary copy of the input with each line assigned a random integer
NUMBERED_LINES=`mktemp` || exit 1
randomIntegers `wc -l < $1` | paste - $1 > $NUMBERED_LINES

# Sort the lines by the first column and then chop that column off
sort -k1,1n $NUMBERED_LINES | cut -f2-

# Remove the temporary file
rm $NUMBERED_LINES
