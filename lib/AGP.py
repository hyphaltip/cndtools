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

class Record:

    def __init__(self, line=None):
        if line is None:
            return

        commentPos = line.find('#')
        if commentPos != -1:
            line = line[:commentPos]

        fields = line.split('\t')
        self.chrom = fields[0]
        self.contigStart, self.contigEnd, self.recNum = map(int, fields[1:4])
        self.type = fields[4]
        if self.type == 'N':
            self.gapLength = int(fields[5])
            self.gapKind = fields[6]
            self.isBridged = (fields[7] == "yes")
        else:
            self.sourceAccession = fields[5]
            self.sourceStart, self.sourceEnd = map(int, fields[6:8])
            self.sourceOrientation = fields[8]

    def isGap(self):
        return self.type == 'N'

    def __str__(self):
        if self.isGap():
            return '\t'.join(map(str, (self.chrom,
                                       self.contigStart,
                                       self.contigEnd,
                                       self.recNum,
                                       self.type,
                                       self.gapLength,
                                       self.gapKind,
                                       (self.isBridged and "yes") or "no")))
        else:
            return '\t'.join(map(str, (self.chrom,
                                       self.contigStart,
                                       self.contigEnd,
                                       self.recNum,
                                       self.type,
                                       self.sourceAccession,
                                       self.sourceStart,
                                       self.sourceEnd,
                                       self.sourceOrientation)))

class Iterator:
    """An iterator over an AGP formatted file.

    Provides an iterator over the records of the AGP file
    """

    def __init__(self, handle):
        self.handle = handle

    def __iter__(self):
        return self

    def next(self):
        while 1:
            line = self.handle.readline()

            # Stop when EOF reached
            if line == "":
                raise StopIteration

            # Skip over comment lines
            if line.startswith('#'):
                continue

            return Record(line[:-1])
