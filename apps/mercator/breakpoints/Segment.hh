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

#ifndef __SEGMENT_HH__
#define __SEGMENT_HH__

#include <string>

#include "types.hh"

class Segment {
private:
	Genome* genome;
	std::string chrom;
	size_t start;
	size_t end;
	char strand;

	SegmentSet* segSet;

	BreakpointSegment* upstreamSegment;
	BreakpointSegment* downstreamSegment;

public:
	Segment(Genome* genome,
			const std::string& chrom,
			size_t start,
			size_t end,
			char strand);

	Genome* getGenome() const;
	std::string getChrom() const;
	size_t getStart() const;
	size_t getEnd() const;
	char getStrand() const;
	size_t getLength() const;

	size_t getUpstreamCoord() const;
	size_t getDownstreamCoord() const;

	void setStart(size_t pos);
	void setEnd(size_t pos);
	
	bool isOnSameChrom(const Segment& s) const;
	bool operator<(const Segment& s) const;

	SegmentSet* getSegSet() const { return segSet; }
	void setSegSet(SegmentSet* segSet) { this->segSet = segSet; }
	
	void setUpstreamSegment(BreakpointSegment* s);
	void setDownstreamSegment(BreakpointSegment* s);
	BreakpointSegment* getUpstreamSegment() const;
	BreakpointSegment* getDownstreamSegment() const;

	friend std::ostream& operator<<(std::ostream& strm, const Segment& s);
};

#endif // __SEGMENT_HH__
