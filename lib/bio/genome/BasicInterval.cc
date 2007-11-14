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

#include <cassert>

#include "bio/genome/BasicInterval.hh"

namespace bio { namespace genome {

	BasicInterval BasicInterval::between(const Interval& i1,
										 const Interval& i2) {
		assert(i1.getChrom() == i2.getChrom());
		assert(i1.getStrand() == i2.getStrand());
		return BasicInterval(i1.getChrom(),
							 std::min(i1.getEnd(), i2.getEnd()),
							 std::max(i1.getStart(), i2.getStart()),
							 i1.getStrand());
	}
	
	BasicInterval::BasicInterval(const Interval& i)
		: chrom(i.getChrom()),
		  start(i.getStart()),
		  end(i.getEnd()),
		  strand(i.getStrand()) {
	}

	BasicInterval::BasicInterval(const std::string& chrom,
								 const Position start,
								 const Position end,
								 const Strand strand) 
		: chrom(chrom), start(start), end(end), strand(strand) {
	}
	
	BasicInterval::BasicInterval(const Coord& start,
								 const Coord& end,
								 const Strand strand)
		: chrom(start.chrom), start(start.pos), end(end.pos), strand(strand) {
		assert(start.chrom == end.chrom);
	}

	BasicInterval::BasicInterval(const Coord& c,
								 const Strand strand)
		: chrom(c.chrom), start(c.pos), end(c.pos), strand(strand) {
	}

	void BasicInterval::setChrom(const std::string& chrom) {
		this->chrom = chrom;
	}
	
	void BasicInterval::setStart(const Position start) {
		this->start = start;
	}
	
	void BasicInterval::setEnd(const Position end) {
		this->end = end;
	}

	void BasicInterval::setStrand(const Strand strand) {
		this->strand = strand;
	}

	std::string BasicInterval::getChrom() const { return chrom; }
	Position BasicInterval::getStart() const { return start; }
	Position BasicInterval::getEnd() const { return end; }
	Strand BasicInterval::getStrand() const { return strand; }
	
	BasicInterval BasicInterval::operator&(const Interval& i) const {
		BasicInterval b = *this;
		b &= i;
		return b;
	}
	
	BasicInterval BasicInterval::operator|(const Interval& i) const {
		BasicInterval b = *this;
		b |= i;
		return b;
	}

	BasicInterval& BasicInterval::operator&=(const Interval& i) {
		assert(chrom == i.getChrom());
		setStart(std::max(getStart(), i.getStart()));
		setEnd(std::min(getEnd(), i.getEnd()));
		return *this;
	}
	
	BasicInterval& BasicInterval::operator|=(const Interval& i) {
		assert(chrom == i.getChrom());
		setStart(std::min(getStart(), i.getStart()));
		setEnd(std::max(getEnd(), i.getEnd()));
		return *this;
	}
	
} }
