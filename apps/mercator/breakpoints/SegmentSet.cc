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

#include "SegmentSet.hh"
#include "Segment.hh"
#include "AlignedSegments.hh"
#include "BreakpointSegment.hh"
#include "Genome.hh"

SegmentSet::SegmentSet() {
}

void SegmentSet::addSegment(Segment* s) {
	segments.push_back(s);
	s->setSegSet(this);
}

void SegmentSet::setNum(size_t num) {
	this->num = num;
}

size_t SegmentSet::getNum() const {
	return num;
}

void SegmentSet::makeAlignedSegments(std::vector<AlignedSegments*>& alignments) {
	typedef std::vector<Segment*>::iterator SegmentIt;
	for (SegmentIt s1It = segments.begin(); s1It != segments.end(); ++s1It) {
		Segment* s1 = *s1It;
		for (SegmentIt s2It = s1It + 1; s2It != segments.end(); ++s2It) {
			Segment* s2 = *s2It;
			// Connect upstream breakpoint segments
			alignments.push_back(new AlignedSegments(s1->getUpstreamSegment(),
													 s2->getUpstreamSegment(),
													 s1->getStrand() == '-',
													 s2->getStrand() == '-'));
			// Connect downstream breakpoint segments
			alignments.push_back(new AlignedSegments(s1->getDownstreamSegment(),
													 s2->getDownstreamSegment(),
													 s1->getStrand() == '+',
													 s2->getStrand() == '+'));
		}
	}
}

void SegmentSet::writeAlignedSegments(std::ostream& stream) {
	typedef std::vector<Segment*>::iterator SegmentIt;

	stream << getNum() << 'u';
	for (SegmentIt it = segments.begin(); it != segments.end(); ++it) {
		Segment* s = *it;

		size_t start = s->getUpstreamSegment()->getBreakpoint();
		size_t end = s->getUpstreamCoord();
		if (s->getStrand() == '+') { std::swap(start, end); }

		stream << '\t' << s->getGenome()->getName()
			   << '\t' << s->getChrom()
			   << '\t' << start
			   << '\t' << end
			   << '\t' << s->getStrand();
	}
	stream << '\n';

	stream << getNum() << 'd';
	for (SegmentIt it = segments.begin(); it != segments.end(); ++it) {
		Segment* s = *it;

		size_t start = s->getDownstreamCoord();
		size_t end = s->getDownstreamSegment()->getBreakpoint();
		if (s->getStrand() == '+') { std::swap(start, end); }

		stream << '\t' << s->getGenome()->getName()
			   << '\t' << s->getChrom()
			   << '\t' << start
			   << '\t' << end
			   << '\t' << s->getStrand();
	}
	stream << '\n';
}

std::ostream& operator<<(std::ostream& strm, const SegmentSet& segmentSet) {
	strm << segmentSet.num;
	typedef std::vector<Segment*>::const_iterator SegmentIt;
	for (SegmentIt it = segmentSet.segments.begin();
		 it != segmentSet.segments.end(); ++it) {
		Segment* s = *it;
		strm << '\t' << s->getGenome()->getName()
			 << '\t' << s->getChrom()
			 << '\t' << s->getStart()
			 << '\t' << s->getEnd()
			 << '\t' << s->getStrand();
	}
	strm << '\n';
	return strm;
}
