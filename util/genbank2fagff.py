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
from Bio import GenBank
from optparse import OptionParser
import GFF
import FASTA

usage = "usage: %prog [options] faFile gffFile < genbankInput"
optparser = OptionParser(usage)
optparser.add_option("--source", dest="source", default="genbank",
                     help="name of the source to use in the gff file")
(options, args) = optparser.parse_args()
if len(args) != 2:
   optparser.error("Invalid number of arguments")

faFile, gffFile = map(file, args, ['w'] * len(args))

parser = GenBank.FeatureParser()
iterator = GenBank.Iterator(sys.stdin, parser)

while 1:
   sys.stderr.write("Parsing record...\n")
   gbrec = iterator.next()

   if gbrec is None:
      break

   seqname = None
   gffRecs = []
   
   for f in gbrec.features:
      if f.type == "source":
         if "chromosome" in f.qualifiers:
            seqname = "chr" + f.qualifiers["chromosome"][0]
         else:
            seqname = gbrec.locus
            
      elif f.type == "CDS":
         if len(f.sub_features) == 0:
            cdss = [(f.location.start.position, f.location.end.position)]
         else:
            cdss = [(sf.location.start.position, sf.location.end.position)
                    for sf in f.sub_features]

         if f.strand == -1:
            cdss.reverse()

         strand = (f.strand == 1 and '+') or '-'
         frame = int(f.qualifiers["codon_start"][0]) - 1

         


         for qual in ("gene", "locus_tag", "protein_id"):
            if qual in f.qualifiers:
               gene_id = f.qualifiers[qual]

         if "gene" in f.qualifiers:
            gene_id = f.qualifiers["gene"][0]
         gene_id = f.qualifiers["gene"][0]


         transcript_id = gene_id
         
         for cds in cdss:
            gffRecs.append(GFF.Record(seqname=seqname,
                                      source=options.source,
                                      feature="CDS",
                                      start=cds[0] + 1,
                                      end=cds[1],
                                      score=None,
                                      strand=strand,
                                      frame=frame,
                                      attributes={"gene_id": gene_id,
                                                  "transcript_id": transcript_id}))
            frame = (cds[1] - cds[0] - frame) % 3
      elif f.type == "mRNA":
         if len(f.sub_features) == 0:
            exons = [(f.location.start.position, f.location.end.position)]
         else:
            exons = [(sf.location.start.position, sf.location.end.position)
                     for sf in f.sub_features]

         if f.strand == -1:
            exons.reverse()

         strand = (f.strand == 1 and '+') or '-'
         gene_id = f.qualifiers["gene"][0]
         transcript_id = gene_id
         
         for exon in exons:
            gffRecs.append(GFF.Record(seqname=seqname,
                                      source=options.source,
                                      feature="exon",
                                      start=exon[0] + 1,
                                      end=exon[1],
                                      score=None,
                                      strand=strand,
                                      frame=None,
                                      attributes={"gene_id": gene_id,
                                                  "transcript_id": transcript_id}))

   faRec = FASTA.Record()
   faRec.sequence = gbrec.seq.tostring()
   faRec.title = seqname

   print >> faFile, faRec

   for gffRec in gffRecs:
      print >> gffFile, gffRec
