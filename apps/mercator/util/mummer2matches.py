#!/usr/bin/env python

import sys

ref_strand = "+"

for line in sys.stdin:
    if line.startswith(">"):
	tokens = line[1:].split()
	assert(1 <= len(tokens) <= 2)
	if len(tokens) == 2:
	    assert(tokens[1] == "Reverse")
	    query_strand = "-"
	else:
	    query_strand = "+"
	query = tokens[0]
    else:
	ref, ref_start, query_start, length = line.split()
	print '\t'.join((ref, ref_strand, ref_start,
			 query, query_strand, query_start,
			 length))
