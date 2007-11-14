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

usage = "usage: %prog [options] genomeName1 genomeName2 ... genomeNameK"
optparser = OptionParser(usage)
optparser.add_option("--genome-dir",
                     dest="genomeDir",
                     help="Path to directory contaning genome sequences " +
                     "in SDB format",
                     default=".", metavar="PATH")
optparser.add_option("--gff-dir", dest="gffDir",
                     help="Path to directory containing gene annotations " +
                     "for each genome in GFF format",
                     default=".", metavar="PATH")
optparser.add_option("--out-dir", dest="outDir",
                     help="Path to directory to output mercator input files",
                     default=".", metavar="PATH")
(options, args) = optparser.parse_args()

genomes = args

if len(genomes) < 2:
    print >>sys.stderr, "Error: At least two genomes must be specified"
    sys.exit(1)

# Check to make sure that a gff and sdb file are present for all genomes
for genome in genomes:
    gffFilename = os.path.join(options.gffDir, "%s.gff" % genome)
    sdbFilename = os.path.join(options.genomeDir, "%s.sdb" % genome)
    
    if not os.path.exists(gffFilename):
        print >>sys.stderr, "Error: GFF file for %s does not exist: %s" % \
              (genome, gffFilename)
        sys.exit(1)
    if not os.path.exists(sdbFilename):
        print >>sys.stderr, "Error: SDB file for %s does not exist: %s" % \
              (genome, sdbFilename)
        sys.exit(1)

# Make anchor, chrom, and protein sequence files
for genome in genomes:
    sdbFilename = os.path.join(options.genomeDir, "%s.sdb" % genome)    
    gffFilename = os.path.join(options.gffDir, "%s.gff" % genome)

    # Make chrom file
    sys.stderr.write("Making chromosome file for %s..." % genome)
    chromFilename = os.path.join(options.outDir, "%s.chroms" % genome)
    os.system("sdbList -l %s > %s" % (sdbFilename, chromFilename))
    sys.stderr.write("done\n")

    # Make anchor file
    sys.stderr.write("Making anchors for %s..." % genome)
    anchorFilename = os.path.join(options.outDir, "%s.anchors" % genome)
    os.system("cat %s | gffRemoveOverlaps -fCDS -l | gff2anchors > %s" % 
	      (gffFilename, anchorFilename))
    sys.stderr.write("done\n")

    # Make protein file
    sys.stderr.write("Extracting protein sequences for anchors...")
    proteinFilename = os.path.join(options.outDir, "%s.proteins.fa" % genome)
    os.system("anchors2fa %s < %s > %s" %
              (sdbFilename, anchorFilename, proteinFilename))
    sys.stderr.write("done\n")

# BLAT protein anchors pairwise
sys.stderr.write("BLATing anchors pairwise...\n")
genomePairs = [(g1, g2) for g1 in genomes for g2 in genomes if g1 < g2]
for (genome1, genome2) in genomePairs:
    sys.stderr.write("%s-%s\n" % (genome1, genome2))
    protein1 = os.path.join(options.outDir, "%s.proteins.fa" % genome1)
    protein2 = os.path.join(options.outDir, "%s.proteins.fa" % genome2)
    blatOutput = os.path.join(options.outDir, "%s-%s.blat" % (genome1, genome2))
    os.system("blat -prot -out=blast8 %s %s %s" %
              (protein2, protein1, blatOutput))
    hitOutput = os.path.join(options.outDir, "%s-%s.hits" % (genome1, genome2))
    os.system("blat2hits < %s > %s" % (blatOutput, hitOutput))
