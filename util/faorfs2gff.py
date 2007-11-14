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
import os
from optparse import OptionParser
import re
import FASTA
import GFF
from Bio.SeqUtils import antiparallel

usage = "usage: %prog dbfile < faorfInput"
optparser = OptionParser(usage)
optparser.add_option("-s", "--source", dest="source", default="source")
(options, args) = optparser.parse_args()
if len(args) != 1:
    optparser.error("Invalid number of arguments")
dbFilename = args[0]

#Example title line:
#ORFN:Sbay_Contig560.7 YLR203C, Contig c560 8347-9657 reverse complement
titlePat = re.compile(r'[^:]*:(\S*) ([^,]*), ' +
                      r'Contig (\S*) (\d+)-(\d+)( reverse complement)?')

def oppositeStrand(strand):
    return (strand == '+' and '-') or '+'

class Genome:
    def __init__(self, dbFilename):
        strm = os.popen("db2fa -c " + dbFilename)
        self.seqs = dict([(r.title, r.sequence) for r in FASTA.Iterator(strm)])
    def getSeq(self, contig, start, end, strand='+'):
        if strand == '+':
            return self.seqs[contig][start:end]
        else:
            return antiparallel(self.seqs[contig][start:end])

g = Genome(dbFilename)

for rec in FASTA.Iterator(sys.stdin):
    m = titlePat.match(rec.title)
    if m is None:
        sys.stderr.write("Could not parse line:\n" + rec.title + '\n')
        continue
    (orf, id, contig, start, end, rc) = m.groups()
    start, end = int(start), int(end)
    strand = (rc is None and '+') or '-'
    opposite = oppositeStrand(strand)
    seq = g.getSeq(contig, start - 1, end, strand)
    rcseq = g.getSeq(contig, start - 1, end, opposite)
    if (rec.sequence == seq):
        pass
    elif (rec.sequence == rcseq):
        strand = oppositeStrand(strand)
        sys.stderr.write("Opposite: " + rec.title + '\n')
    else:
        sys.stderr.write("No Match: " + rec.title + '\n')
        continue

    gffRec = GFF.Record(seqname=contig,
                        source=options.source,
                        feature="CDS",
                        start=start,
                        end=end,
                        score=None,
                        strand=strand,
                        frame=0,
                        attributes={"gene_id": id or orf})
    print gffRec
