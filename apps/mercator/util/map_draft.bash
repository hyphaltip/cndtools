#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: $0 draft_genome target_genome draft_agpfile align_dir"
    exit 1
fi

QUERY=$1
TARGET=$2
AGP_FILE=$3
ALIGN_DIR=$4

# Temporary files for the target and query genome intervals (in GFF format)
TARGET_GFF=`mktemp` || exit 1
QUERY_GFF=`mktemp` || exit 1

# Map the contigs (in assembled coordinates) to the target genome
awk '$5 == "D" {OFS="\t"; print $1,".",$6,$2,$3,".",$9,".","";}' < $AGP_FILE | \
    gffMap $ALIGN_DIR $QUERY $TARGET > $TARGET_GFF

# Map the target genome intervals back to the draft genome contigs
gffMap --output-unmapped $ALIGN_DIR $TARGET $QUERY < $TARGET_GFF | \
    gffUntransform --ignore $AGP_FILE > $QUERY_GFF

# Make sure the number of lines is the same in both files
QUERY_LINES=`wc -l $QUERY_GFF | awk '{print $1}'`
TARGET_LINES=`wc -l $TARGET_GFF | awk '{print $1}'`
if [ $QUERY_LINES -ne $TARGET_LINES ]; then
    echo "Error: Query and target have a different number of mapped lines"
    exit 1
fi

# Output the two sets of intervals
paste $QUERY_GFF $TARGET_GFF | grep -v UNMAPPED | cut -f1,4,5,7,10,13,14,16

# Remove the temporary files
rm $TARGET_GFF $QUERY_GFF
