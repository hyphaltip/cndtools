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

#ifndef __BREAKPOINT_SEGMENT_HH__
#define __BREAKPOINT_SEGMENT_HH__

#include <string>
#include <vector>

#include "types.hh"

class BreakpointSegment {
private:
	static size_t resolution;
	
	Segment* upstreamSegment;
	Segment* downstreamSegment;
	size_t num;
	bool changed;
	
	std::vector<Score> scores;
	size_t bestIndex;

	size_t lowerBound;
	size_t upperBound;
	
public:
	static void setResolution(size_t resolution);
	static size_t getResolution();

	
	BreakpointSegment(Segment* upstreamSegment,
					  Segment* downstreamSegment);

	Genome* getGenome() const;
	std::string getChrom() const;
	size_t getStart() const;
	size_t getEnd() const;

	bool isColinearWith(const BreakpointSegment& other) const;
	Segment* getUpstreamSegment() const { return upstreamSegment; }
	Segment* getDownstreamSegment() const { return downstreamSegment; }
	
	// Returns true if the bounds of this breakpoint segment changed
	// during the last shrink operation
	bool isChanged() const;
	
	std::string getSeq(char strand = '+') const;

	void setNum(size_t num);
	size_t getNum() const;

	size_t getLength() const;
	size_t getBoundedLength() const;

	void initIndices(size_t maxSegmentPositions);
	
	size_t getPosition(unsigned long long i) const;

	Score getScore(size_t i) const;
	void setScore(size_t i, Score score);

	size_t getNumIndices() const;
	size_t getBestIndex() const;
	void setBestIndex(size_t i);

	void initBounds();
// 	void updateBounds();
	void updateBounds(size_t minBoundedLength);
	
	void resetScores();
	
	size_t getBreakpoint() const;
	size_t getRelBreakpoint() const;
	void setRelBreakpoint(size_t bp);
	void setBestBreakpoint();
	
	// Split this segment at the breakpoint
	void split();
	
	friend std::ostream& operator<<(std::ostream& strm, const BreakpointSegment& seg);
};

#endif // __BREAKPOINT_SEGMENT_HH__
