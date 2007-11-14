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

from __future__ import generators

def strToInt(field):
    if field == '.':
        return 0
    else:
        return int(field)

def intToStr(field):
    if field == 0:
        return '.'
    else:
        return str(field)

class Record:
    def __init__(self,
                 seqname="",
                 features=None):
        self.seqname = seqname
        if features is None:
            self.features = []
        else:
            self.features = features

    def __str__(self):
        return '\n'.join([">%s" % self.seqname] + map(str, self.features))

class Feature:
    def __init__(self,
                 line=None,
                 label="",
                 start=0,
                 end=0,
                 group=""):
        if line is None:
            self.label = label
            self.start = begin
            self.end = end
            self.group = group
        else:
            fields = line[:-1].split('\t')
            self.label = fields[0]
            self.start = int(fields[1])
            self.end = int(fields[2])
            self.group = fields[3]

    def __str__(self):
        return '\t'.join((self.label,
                          str(self.start),
                          str(self.end),
                          self.group))

class LongFeature(Feature):
    def __init__(self,
                 line=None,
                 label="",
                 start=0,
                 end=0,
                 strand=".",
                 score=0.0,
                 fivePrimeOverhang=0,
                 threePrimeOverhang=0,
                 frame=0,
                 group=""):
        
        if line is None:
            self.label = label
            self.start = start
            self.end = end
            self.strand = strand
            self.score = score
            self.fivePrimeOverhang = fivePrimeOverhang
            self.threePrimeOverhang = threePrimeOverhang
            self.frame = frame
            self.group = group
        else:
            fields = line[:-1].split('\t')
            self.label = fields[0]
            self.start = int(fields[1])
            self.end = int(fields[2])
            self.strand = fields[3]
            self.score = float(fields[4])
            self.fivePrimeOverhang = strToInt(fields[5])
            self.threePrimeOverhang = strToInt(fields[6])
            self.frame = strToInt(fields[7])
            self.group = fields[8]

    def __str__(self):
        return '\t'.join((self.label,
                          str(self.start),
                          str(self.end),
                          self.strand,
                          str(self.score),
                          intToStr(self.fivePrimeOverhang),
                          intToStr(self.threePrimeOverhang),
                          intToStr(self.frame),
                          self.group))

def Iterator(strm):
    while 1:
        line = strm.readline()
        if not line:
            return
        elif line.startswith('>'):
            rec = Record(line[1:-1])
            break
    
    while 1:
        line = strm.readline()
        if not line:
            yield rec
            return
        elif line.startswith('>'):
            yield rec
            rec = Record(line[1:-1])
        elif line.count('\t') == 3:
            rec.features.append(Feature(line=line))
        else:
            rec.features.append(LongFeature(line=line))
