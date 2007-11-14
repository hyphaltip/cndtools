/* Copyright (c) 2006
   Colin Dewey (University of Wisconsin-Madison)
   cdewey@biostat.wisc.edu
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "Segment.hh"
#include "Genome.hh"

Segment::Segment(Genome* genome,
				 const std::string& chrom,
				 size_t start,
				 size_t end,
				 char strand)
	: genome(genome), chrom(chrom), start(start), end(end), strand(strand)
{}

Genome* Segment::getGenome() const { return genome; }
std::string Segment::getChrom() const { return chrom; }
size_t Segment::getStart() const { return start; }
size_t Segment::getEnd() const { return end; }
char Segment::getStrand() const { return strand; }

size_t Segment::getLength() const { return end - start; }

bool Segment::isOnSameChrom(const Segment& s) const {
	return genome->getNum() == s.genome->getNum() && chrom == s.chrom;
}

bool Segment::operator<(const Segment& s) const {
	return genome->getNum() < s.genome->getNum()
		|| (genome->getNum() == s.genome->getNum() &&
			(chrom < s.chrom || (chrom == s.chrom && start < s.start)));
}			

size_t Segment::getUpstreamCoord() const {
	return (strand == '+' ? start : end);
}

size_t Segment::getDownstreamCoord() const {
	return (strand == '+' ? end : start);
}

void Segment::setUpstreamSegment(BreakpointSegment* s) {
	upstreamSegment = s;
}

void Segment::setDownstreamSegment(BreakpointSegment* s) {
	downstreamSegment = s;
}

BreakpointSegment* Segment::getUpstreamSegment() const {
	return upstreamSegment;
}

BreakpointSegment* Segment::getDownstreamSegment() const {
	return downstreamSegment;
}

void Segment::setStart(size_t pos) {
	start = pos;
}

void Segment::setEnd(size_t pos) {
	end = pos;
}

std::ostream& operator<<(std::ostream& strm, const Segment& s) {
	strm << '(' << s.getGenome()->getName()
		 << ':' << s.getChrom()
		 << ':' << s.getStart()
		 << '-' << s.getEnd()
		 << ':' << s.getStrand();
	return strm;
}
