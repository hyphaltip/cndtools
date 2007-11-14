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


from optparse import OptionParser
import sys

import ZFF
import GFF

usage = "usage: %prog < zffInput > gtfOutput"
optparser = OptionParser(usage)
optparser.add_option("--source", dest="source", default="source")
optparser.add_option("--seglength", type="int", dest="segLength",
                     help="sequences were broken into segments of size LEN",
                     metavar="LEN", default=0)
(options, args) = optparser.parse_args()
if len(args) != 0:
    optparser.error("Invalid number of arguments")

source = options.source

exonTypes = {"Einit": "initial",
             "Exon": "internal",
             "Eterm": "terminal",
             "Esngl": "single"}

for rec in ZFF.Iterator(sys.stdin):
    if options.segLength != 0:
        splitIndex = rec.seqname.rfind(".")
        segNum = int(rec.seqname[splitIndex + 1:])
        offset = segNum * options.segLength
        rec.seqname = rec.seqname[:splitIndex]
    else:
        offset = 0
    
    for feature in rec.features:
        feature.start += offset
        feature.end += offset
        
        if feature.label not in exonTypes:
            print GFF.Record(seqname=rec.seqname,
                             source=source,
                             feature=feature.label,
                             start=feature.start,
                             end=feature.end,
                             score=feature.score,
                             strand=feature.strand,
                             frame=None,
                             attributes={"group": [feature.group]})
        else:
            exontype = exonTypes[feature.label]
            gene_id = feature.group
            transcript_id = "%s.1" % gene_id

            frame = feature.fivePrimeOverhang

            if exontype in ("initial", "single"):
                if feature.strand == '+':
                    startCodonStart = feature.start
                    startCodonEnd = feature.start + 2
                else:
                    startCodonStart = feature.end - 2
                    startCodonEnd = feature.end
                
                print GFF.Record(seqname=rec.seqname,
                                 source=source,
                                 feature="start_codon",
                                 start=startCodonStart,
                                 end=startCodonEnd,
                                 score=None,
                                 strand=feature.strand,
                                 frame=None,
                                 attributes={"gene_id": [gene_id],
                                             "transcript_id": [transcript_id]})
            
            print GFF.Record(seqname=rec.seqname,
                             source=source,
                             feature='CDS',
                             start=feature.start,
                             end=feature.end,
                             score=feature.score,
                             strand=feature.strand,
                             frame=frame,
                             attributes={"gene_id": [gene_id],
                                         "transcript_id": [transcript_id],
                                         "exontype": [exontype]})

            print GFF.Record(seqname=rec.seqname,
                             source=source,
                             feature='exon',
                             start=feature.start,
                             end=feature.end,
                             score=feature.score,
                             strand=feature.strand,
                             frame=None,
                             attributes={"gene_id": [gene_id],
                                         "transcript_id": [transcript_id],
                                         "exontype": [exontype]})

            if exontype in ("terminal", "single"):
                if feature.strand == '+':
                    stopCodonStart = feature.end -2
                    stopCodonEnd = feature.end
                else:
                    stopCodonStart = feature.start
                    stopCodonEnd = feature.start + 2
                
                print GFF.Record(seqname=rec.seqname,
                                 source=source,
                                 feature="stop_codon",
                                 start=stopCodonStart,
                                 end=stopCodonEnd,
                                 score=None,
                                 strand=feature.strand,
                                 frame=None,
                                 attributes={"gene_id": [gene_id],
                                             "transcript_id": [transcript_id]})
