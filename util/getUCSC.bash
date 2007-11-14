#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 UCSCGenome table [test=0]"
    exit
fi

DB=$1
TABLE=${DB}.$2
if [ ${3:-"0"} -eq "0" ]; then
    URL="http://genome.ucsc.edu/cgi-bin/hgText"
else
    URL="http://genome-test.cse.ucsc.edu/cgi-bin/hgText"
fi

OUTPUT_TYPE="Tab-separated%2C%20All%20fields"

EXTRA_ARGS="&position=genome&tbPosOrKeys=pos&outputType=${OUTPUT_TYPE}&phase=Get results"

wget -q -O - "${URL}?db=${DB}&table=${TABLE}${EXTRA_ARGS}" | perl -p -e "s/$1_$2/$2/"
