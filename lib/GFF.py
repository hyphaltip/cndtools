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
import re

class Record:
    """A record from a GFF file.

    Fields:
       seqname
       source
       feature
       start
       end
       score
       strand
       frame
       attributes
       comments
    """

    fieldNames = ["seqname", "source", "feature", "start", "end", "score",
                  "strand", "frame", "attributes", "comments"]
    
    def __init__(self, gffLine=None, **keywords):
        if gffLine is None:
            self.__dict__.update(dict(zip(Record.fieldNames,
                                          [None] * len(Record.fieldNames))))
            self.__dict__.update(keywords)
            if self.attributes is None:
                self.attributes = {}
            return

        fields = gffLine[:-1].split('\t', 8)
        
        # Check for proper field structure
        if len(fields) not in (8, 9):
            raise Exception, "Invalid GFF line: %s" % gffLine
        
        self.seqname = fields[0]
        self.source = fields[1]
        self.feature = fields[2]
        self.start = long(fields[3])
        self.end = long(fields[4])
        
        if fields[5] == '.':
            self.score = None
        else:
            self.score = float(fields[5])
                
        if fields[6] == '.':
            self.strand = None
        else:
            self.strand = fields[6]
            
        if fields[7] == '.':
            self.frame = None
        else:
            self.frame = int(fields[7])

        self.attributes = {}

        if fields[8]:
            self.parseAttributes(fields[8])

    def parseAttributes(self, s):
        currentTag = None
        for token in AttributeIterator(s):
            if currentTag is None:
                if isinstance(token, IdentifierToken):
                    currentTag = token.value
                    self.attributes[currentTag] = []
                else:
                    raise Exception, "Invalid GFF attributes: " + s
            elif isinstance(token, SeparatorToken):
                currentTag = None
            elif isinstance(token, CommentToken):
                return
            elif isinstance(token, IdentifierToken):
                self.attributes[currentTag].append(token.value)
            elif isinstance(token, ValueToken):
                self.attributes[currentTag].append(token.value)
            else:
                raise Exception, "Invalid GFF attributes: " + s

    def copy(self):
        return Record(seqname=self.seqname,
                      source=self.source,
                      feature=self.feature,
                      start=self.start,
                      end=self.end,
                      score=self.score,
                      strand=self.strand,
                      frame=self.frame,
                      attributes=self.attributes.copy())

    def __repr__(self):
        return repr(self.__dict__)

    def __str__(self):
        fields = [(f is None and '.') or str(f) for f in (self.seqname,
                                                          self.source,
                                                          self.feature,
                                                          self.start,
                                                          self.end,
                                                          self.score,
                                                          self.strand,
                                                          self.frame)]
        attrs = ' '.join([' '.join([tag] + map(quote, values)) + ";"
                          for tag, values in self.attributes.items()])
        fields.append(attrs)
        return '\t'.join(fields)

def quote(s):
    return '"%s"' % str(s)

class IdentifierToken:
    pass
class ValueToken:
    pass
class CommentToken:
    pass
class SeparatorToken:
    pass
class UnknownToken:
    pass

class AttributeIterator:
    identifierPat = re.compile(r'\s*([A-Za-z][A-Za-z0-9_]*)')
    freeTextPat = re.compile(r'\s*"(([^"]|(\\"))*)(?<!\\)"')
    valuePat = re.compile(r'\s*([^;# \t\n\r\f\v]+)')
    sepPat = re.compile(r'\s*(;)')
    commentPat = re.compile(r'\s*#(.*)$')

    pats = (identifierPat, freeTextPat, valuePat, sepPat, commentPat)
    tokenClasses = (IdentifierToken, ValueToken, ValueToken,
                    SeparatorToken, CommentToken)

    def __init__(self, s):
        self.s = s.rstrip()
        self.pos = 0
        
    def __iter__(self):
        return self
    
    def next(self):
        if self.pos >= len(self.s):
            raise StopIteration
        
        for (pat, tclass) in zip(AttributeIterator.pats,
                                 AttributeIterator.tokenClasses):
            match = pat.match(self.s, self.pos)
            if match is not None:
                self.pos = match.end(0)
                t = tclass()
                t.value = match.group(1)
                return t
        else:
            return UnknownToken()

class Iterator:
    """An iterator over a GFF formatted file.

    Provides an iterator over the records of the GFF file
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

            return Record(line)
