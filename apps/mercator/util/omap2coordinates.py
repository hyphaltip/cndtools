#!/usr/bin/env python

import sys
import os
import AGP
from optparse import OptionParser

usage = "usage: %prog [options] < map"
optparser = OptionParser(usage)
optparser.add_option("--map-dir", dest="map_dir",
                     help="Directory containing map output files (default: .)",
                     default=".", metavar="DIR")
(options, args) = optparser.parse_args()
if len(args) != 0:
    optparser.error("incorrect number of arguments")

def read_agp_file(stream):
    recs = {}
    for rec in AGP.Iterator(stream):
	recs.setdefault(rec.chrom, []).append(rec)
    return recs

def open_agp_files(genomes, mapDir):
    agp_files = {}
    for g in genomes:
	agp_filename = os.path.join(mapDir, g + ".agp")
	if os.path.exists(agp_filename):
	    agp_files[g] = read_agp_file(file(agp_filename))
    return agp_files

def make_interval(fields):
    if fields[0] == "NA":
	return None
    else:
	return [fields[0], int(fields[1]), int(fields[2]), fields[3]]

def split_intervals(fields):
    return [make_interval(fields[i: i + 4]) for i in xrange(0, len(fields), 4)]

def transform_interval(interval, agp_recs):
    transformed = []
    for rec in agp_recs:
	if interval[1] >= rec.contigEnd:
	    continue
	elif interval[2] < rec.contigStart:
	    break
	elif rec.isGap():
	    continue
	else:
	    rec_len = rec.contigEnd - (rec.contigStart - 1)
	    start_offset = max(0, interval[1] - (rec.contigStart - 1))
	    end_offset = min(rec_len, interval[2] - (rec.contigStart - 1))
	    if rec.sourceOrientation == interval[3]:
		orientation = '+'
	    else:
		orientation = '-'
	    
	    t = [rec.sourceAccession, 0, 0, orientation]
	    if rec.sourceOrientation == '+':
		t[1] = (rec.sourceStart - 1) + start_offset
		t[2] = (rec.sourceStart - 1) + end_offset
	    else:
		t[1] = rec.sourceEnd - end_offset
		t[2] = rec.sourceEnd - start_offset
	    transformed.append(t)

    if interval[3] == '-':
	transformed.reverse()
    
    return transformed

def tab_join(fields):
    return '\t'.join(map(str, fields))

genomes = file(os.path.join(options.map_dir, "genomes")).read().split()
draft_agps = open_agp_files(genomes, options.map_dir)

for line in sys.stdin:
    fields = line.split()
    num = fields[0]
    intervals = split_intervals(fields[1:])
    for g, i in zip(genomes, intervals):
	if i is None: 
	    continue
	elif g in draft_agps:
	    transformed = transform_interval(i, draft_agps[g][i[0]])
	    print tab_join([num, g] + [tab_join(t) for t in transformed])
	else:
	    print tab_join([num, g] + i)
