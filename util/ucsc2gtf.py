#!/usr/bin/env python

# Copyright (c) 2006
# Colin Dewey (University of Wisconsin-Madison)
# cdewey@biostat.wisc.edu
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


import sys
import GFF

if len(sys.argv) != 2:
    sys.stderr.write("Usage: ucsc2gtf sourceName < ucscInput\n")
    sys.exit(1)

source = sys.argv[1]

# Dictionary to keep track of transcript numbers
gene_ids = {}

for line in sys.stdin:
    fields = line[:-1].split('\t')

    if len(fields) == 16:
        fields = fields[1:]
    
    (name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, \
     exonCount, exonStarts, exonEnds) = fields[0:10]

    (txStart, txEnd, cdsStart, cdsEnd, exonCount) = \
              map(long, (txStart, txEnd, cdsStart, cdsEnd, exonCount))

    transcript_id = "%s.%d" % (name, gene_ids.setdefault(name, 1))
    gene_ids[name] += 1

    rec = GFF.Record(seqname=chrom,
                     source=source,
                     feature=None,
                     start=None,
                     end=None,
                     score=None,
                     strand=strand,
                     frame=None,
                     attributes={"gene_id": [name],
                                 "transcript_id": [transcript_id]})

    features = []

    rec.feature = "start_codon"
    if strand == '+':
        rec.start = cdsStart + 1
        rec.end = cdsStart + 3
    else:
        rec.start = cdsEnd - 2
        rec.end = cdsEnd
    features.append(str(rec))
        

    # Generate gff lines for exons
    exonStarts = map(long, exonStarts.split(',')[:-1])
    exonEnds = map(long, exonEnds.split(',')[:-1])
    if strand == '-':
        exonStarts.reverse()
        exonEnds.reverse()
        
    frame = 0
    for (exonNum, exonStart, exonEnd) in zip(xrange(exonCount),
                                             exonStarts,
                                             exonEnds):
        if exonNum == 0:
            rec.attributes["exontype"] = ["initial"]
        elif exonNum == (exonCount - 1):
            rec.attributes["exontype"] = ["terminal"]
        else:
            rec.attributes["exontype"] = ["internal"]

        rec.feature = "exon"
        rec.start = exonStart + 1
        rec.end = exonEnd
        rec.frame = None
        features.append(str(rec))

        exonCodingStart = max(exonStart, cdsStart)
        exonCodingEnd = min(exonEnd, cdsEnd)
        if (exonCodingStart < exonCodingEnd):
            rec.feature = "CDS"
            rec.start = exonCodingStart + 1
            rec.end = exonCodingEnd
            rec.frame = (3 - frame) % 3
            frame = (frame + exonCodingEnd - exonCodingStart) % 3
            features.append(str(rec))                    

    rec.feature = "stop_codon"
    rec.frame = None
    del rec.attributes["exontype"]
    if strand == '+':
        rec.start = cdsEnd - 2
        rec.end = cdsEnd
    else:
        rec.start = cdsStart + 1
        rec.end = cdsStart + 3
    features.append(str(rec))        

    # Write GFF lines for this gene
    sys.stdout.write('\n'.join(features))
    sys.stdout.write('\n')
