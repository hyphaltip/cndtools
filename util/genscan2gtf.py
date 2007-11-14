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
import AGP
from optparse import OptionParser

usage = "usage: %prog [options] < genscanInput > gtfOutput"
parser = OptionParser(usage)
parser.add_option("-s", "--seglength", type="int", dest="segLength",
                  help="sequences were broken into segments of size LEN",
                  metavar="LEN", default=0)
parser.add_option("-a", "--agpfile", dest="agpFilename",
                  help="genscan was run on contigs of an assembly.  Convert coordinates based on the agp file FILENAME",
                  metavar="FILENAME")
(options, args) = parser.parse_args()
if len(args) != 0:
    parser.error("incorrect number of arguments")

class SegmentMapper:
    def __init__(self, segLength=0, agpFilename=None):
        self.segLength = segLength
        self.segMap = {}
        if agpFilename is not None:
            sys.stderr.write("Reading AGP file...")
            for rec in AGP.Iterator(file(agpFilename)):
                if not rec.isGap():
                    assert(rec.sourceStart == 1)
                    self.segMap[rec.sourceAccession] = (rec.chrom,
                                                        rec.contigStart,
                                                        rec.contigEnd,
                                                        rec.sourceOrientation)
            sys.stderr.write("done\n")
    
    def getPlacement(self, segname):
        if self.segLength != 0:
            segNum = int(segname[segname.rfind(".") + 1:])
            offset = segNum * self.segLength
            seqname = segname[:segname.rfind(".")]
        else:
            offset = 0
            seqname = segname

        if seqname in self.segMap:
            placement = self.segMap[seqname]
            if placement[3] == '+':
                return (placement[0], placement[1] + offset, 1)
            else:
                return (placement[0], placement[2] - offset, -1)
        else:
            return (seqname, 1 + offset, 1)

class NotBlankLineException(Exception):
    pass

class EndOfFileException(Exception):
    pass

def consumeLinesUntilStartsWith(handle, lineStart):
    while 1:
        line = handle.readline()
        if not line:
            raise EndOfFileException
        if line.startswith(lineStart):
            return line

def consumeBlankLines(handle, numLines):
    for i in xrange(numLines):
        line = handle.readline()
        if not line:
            raise EndOfFileException
        if line != "\n":
            raise NotBlankLineException, line

exonTypes = {"Init": "initial",
             "Intr": "internal",
             "Term": "terminal",
             "Sngl": "single"}

def consumePredictions(handle, segname, seqname,
                       offset=1, sign=1, source="genscan"):
    recs = []
    suboptimal = 0
    suboptNum = 1
    while 1:
        line = handle.readline()
        if not line:
            raise EndOfFileException
        elif line == "\n":
            continue
        elif line.startswith("Suboptimal"):
            suboptimal = 1
            consumeLinesUntilStartsWith(handle, "-----")
            continue
        elif line.startswith("NO"):
            continue
        elif line.startswith("Predicted"):
            return recs

        vals = line.rstrip().split()

        if vals[1] in exonTypes:
            id, type, segStrand = vals[0:3]
            begin, end = map(int, vals[3:5])
            if vals[11].endswith("+"): # Handle 0.999+ case
                score = float(vals[11][:-1])
            else:
                score = float(vals[11])
            
            exontype = exonTypes[type]

            geneNum, exonNum = id.split('.')

            if sign == 1:
                strand = segStrand
            else:
                strand = (segStrand == '+' and '-') or '+'

            if suboptimal:
                gene_id = "%s.S%04d" % (segname, suboptNum)
                suboptNum += 1
            else:
                gene_id = "%s.%04d" % (segname, int(geneNum))  

            transcript_id = "%s.1" % gene_id

            if segStrand == '+':
                frame = (int(vals[6]) - begin + 1) % 3
            else:
                frame = (begin - int(vals[6])) % 3

            begin -= 1
            end -= 1

            if exontype == "initial":
                if segStrand == '+':
                    startCodonStart = begin
                    startCodonEnd = begin + 2
                else:
                    startCodonEnd = begin
                    startCodonStart = begin - 2

                if sign == -1:
                    startCodonStart, startCodonEnd = startCodonEnd, startCodonStart
                
                recs.append(GFF.Record(seqname=seqname,
                                       source=source,
                                       feature="start_codon",
                                       start=startCodonStart * sign + offset,
                                       end=startCodonEnd * sign + offset,
                                       score=None,
                                       strand=strand,
                                       frame=None,
                                       attributes={"gene_id": [gene_id],
                                                   "transcript_id": [transcript_id]}))

            recs.append(GFF.Record(seqname=seqname,
                                   source=source,
                                   feature="CDS",
                                   start=min(begin * sign, end * sign) + offset,
                                   end=max(begin * sign, end * sign) + offset,
                                   score=score,
                                   strand=strand,
                                   frame=frame,
                                   attributes={"gene_id": [gene_id],
                                               "transcript_id": [transcript_id],
                                               "exontype": [exontype]}))

            recs.append(GFF.Record(seqname=seqname,
                                   source=source,
                                   feature="exon",
                                   start=min(begin * sign, end * sign) + offset,
                                   end=max(begin * sign, end * sign) + offset,
                                   score=score,
                                   strand=strand,
                                   frame=None,
                                   attributes={"gene_id": [gene_id],
                                               "transcript_id": [transcript_id],
                                               "exontype": [exontype]}))
            
            if exontype == "terminal":
                if segStrand == '+':
                    stopCodonStart = end - 2
                    stopCodonEnd = end
                else:
                    stopCodonEnd = end + 2
                    stopCodonStart = end

                if sign == -1:
                    stopCodonStart, stopCodonEnd = stopCodonEnd, stopCodonStart
                
                recs.append(GFF.Record(seqname=seqname,
                                       source=source,
                                       feature="stop_codon",
                                       start=stopCodonStart * sign + offset,
                                       end=stopCodonEnd * sign + offset,
                                       score=None,
                                       strand=strand,
                                       frame=None,
                                       attributes={"gene_id": [gene_id],
                                                   "transcript_id": [transcript_id]}))

genscanInput = sys.stdin

segMapper = SegmentMapper(segLength=options.segLength,
                          agpFilename=options.agpFilename)

while 1:
    try:
        seqLine = consumeLinesUntilStartsWith(genscanInput, "Sequence")
        segname = seqLine[9:seqLine.find(" : ")]

        tableLine = consumeLinesUntilStartsWith(genscanInput, "----- ")
        consumeBlankLines(genscanInput, 1)

        placement = segMapper.getPlacement(segname)
        
        gffRecs = consumePredictions(genscanInput,
                                     segname,
                                     seqname=placement[0],
                                     offset=placement[1],
                                     sign=placement[2])
        for rec in gffRecs:
            print rec
    except EndOfFileException:
        break
    except NotBlankLineException, e:
        print "Error parsing file:", str(e)
        break
