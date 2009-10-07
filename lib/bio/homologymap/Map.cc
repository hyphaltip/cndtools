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

#include "bio/homologymap/Map.hh"
#include "util/io/line/InputStream.hh"
#include <algorithm>

namespace bio { namespace homologymap {

	bool SegmentNumSorter::operator()(const Segment* s1, const Segment* s2) const {
		return s1->num < s2->num;
	}

	SegmentCoordSorter::SegmentCoordSorter(const size_t genomeNum)
		: genomeNum(genomeNum) {
	}

	bool SegmentCoordSorter::operator()(const Segment* s1, const Segment* s2) const {
		return (*s1->intervals[genomeNum]) < (*s2->intervals[genomeNum]);
	}

	Map::Map() :
		segments(),
		sortedSegments() {
	}

	void Map::read(std::istream& strm) {
		// Read all of the segments in
		util::io::line::InputStream lineStream(strm);
		std::string line;
		Segment seg;
		while (lineStream >> line) {
			line >> seg;
			segments.push_back(new Segment(seg));
		}

		// There is nothing to be done if there are no segments
		if (segments.empty()) {
			return;
		}

		// Sort segments by segment number
		std::sort(segments.begin(), segments.end(), SegmentNumSorter());

		// Initialize sorted segment vectors
		size_t numGenomes = segments[0]->intervals.size();
		sortedSegments = std::vector< std::vector<Segment*> >(numGenomes);

		// Fill sorted segment vectors
		for (std::vector<Segment*>::const_iterator pos = segments.begin();
			 pos != segments.end(); ++pos) {
			for (size_t g = 0; g < numGenomes; ++g) {
				if ((*pos)->hasGenome(g)) {
					sortedSegments[g].push_back(*pos);
				}
			}
		}

		// Sort sorted segment vectors
		for (size_t g = 0; g < numGenomes; ++g) {
			std::sort(sortedSegments[g].begin(), sortedSegments[g].end(),
					  SegmentCoordSorter(g));
		}
	}

	void Map::write(std::ostream& strm) const {
		for (std::vector<Segment*>::const_iterator pos = segments.begin();
			 pos != segments.end(); ++pos) {
			strm << **pos;
		}
	}

	size_t Map::getNumSegments() const {
		return segments.size();
	}
	
	size_t Map::getNumSegments(size_t genomeNum) const {
		return sortedSegments[genomeNum].size();
	}

	Segment* Map::getSegment(const size_t num) const {
		return segments[num];
	}
			
	Segment* Map::getSegment(const size_t num, const size_t genomeNum) const {
		return sortedSegments[genomeNum][num];
	}

	Segment* Map::getSegment(const size_t genomeNum,
							 const std::string& chrom,
							 const genome::Position pos) const {
		genome::Coord coord(chrom, pos);
		return getSegment(genomeNum, coord);
	}
	
	Segment* Map::getSegment(const size_t genomeNum,
							 const genome::Coord& coord) const {
		std::vector<Segment*>::const_iterator pos;
		pos = getSegmentIter(genomeNum, coord);
		if (pos == sortedSegments[genomeNum].end() ||
			!((*pos)->intervals[genomeNum]->contains(coord))) {
			return NULL;
		} else {
			return *pos;
		}
	}

	std::vector<Segment*>::const_iterator
	Map::getSegmentIter(const size_t genomeNum,
						const genome::Coord& coord) const {
		// Do a binary search on the sorted segments for this genome
		typedef std::vector<Segment*>::const_iterator Iter;
		Iter begin = sortedSegments[genomeNum].begin();
		Iter end = sortedSegments[genomeNum].end();
		while (begin != end) {
			Iter middle = begin + (end - begin) / 2;
			genome::Interval* midInterval = (*middle)->intervals[genomeNum];
			if (midInterval->contains(coord)) {
				return middle;
			} else if (coord < *midInterval) {
				end = middle;
			} else {
				begin = middle + 1;
			}
		}
		
		return begin;
	}

	void Map::getSegments(const size_t genomeNum,
						  const genome::Interval& genomeInt,
						  std::vector<Segment*>& segs) const {
		genome::Coord startCoord = genomeInt.getStartCoord();
		genome::Coord endCoord = genomeInt.getEndCoord();
		
		std::vector<Segment*>::const_iterator pos;
		for (pos = getSegmentIter(genomeNum, startCoord);
			 (pos != sortedSegments[genomeNum].end() &&
			  (*pos)->intervals[genomeNum]->getStartCoord() < endCoord);
			 ++pos) {
			segs.push_back(*pos);
		}
	}

} }
