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

#ifndef __BIO_HOMOLOGYMAP_MAP_HH__
#define __BIO_HOMOLOGYMAP_MAP_HH__

#include "bio/homologymap/Segment.hh"

namespace bio { namespace homologymap {

    class Map {
	public:
		Map();

		void read(std::istream& strm);
		void write(std::ostream& strm) const;

		Segment* getSegment(const size_t num) const;
		Segment* getSegment(const size_t num, const size_t genomeNum) const;
		Segment* getSegment(const size_t genomeNum,
							const genome::Coord& coord) const;
		Segment* getSegment(const size_t genomeNum,
							const std::string& chrom,
							const genome::Position pos) const;

		void getSegments(const size_t genomeNum,
						 const genome::Interval& genomeInt,
						 std::vector<Segment*>& segs) const;
		
		size_t getNumSegments() const;
		size_t getNumSegments(size_t genomeNum) const;
	private:
		std::vector<Segment*> segments;
		std::vector< std::vector<Segment*> > sortedSegments;

		std::vector<Segment*>::const_iterator
		getSegmentIter(const size_t genomeNum,
					   const genome::Coord& coord) const;

	};

	struct SegmentNumSorter {
		bool operator()(const Segment* s1, const Segment* s2) const;
	};

	struct SegmentCoordSorter {
		size_t genomeNum;
		SegmentCoordSorter(const size_t genomeNum);
		bool operator()(const Segment* s1, const Segment* s2) const;
	};

} }

#endif // __BIO_HOMOLOGYMAP_MAP_HH__
