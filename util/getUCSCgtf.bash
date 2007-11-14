#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 UCSCGenome table [test=0]"
    exit
fi

DB=$1
TABLE=$2
if [ ${3:-"0"} -eq "0" ]; then
    URL="http://genome.ucsc.edu/cgi-bin/hgTables"
else
    URL="http://genome-test.cse.ucsc.edu/cgi-bin/hgTables"
fi

EXTRA_ARGS="&hgta_regionType=genome&hgta_outputType=gff&&hgta_compressType=hgta_compressNone&hgta_doTopSubmit=get+output"
URL="${URL}?db=${DB}&hgta_track=${TABLE}&hgta_table=${TABLE}${EXTRA_ARGS}"

wget -q -O - "${URL}" | perl -p -e "s/$1_$2/$2/"
