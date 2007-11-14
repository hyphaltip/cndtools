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

#include "bio/genome/Interval.hh"

namespace bio { namespace genome {

	Interval::~Interval() {}
	
	bool Interval::contains(const Coord& c) const {
		return c.pos >= getStart()
			and c.pos < getEnd()
			and c.chrom == getChrom();
	}

	bool Interval::contains(const Interval& i) const {
		return i.getStart() >= getStart()
			and i.getEnd() <= getEnd()
			and i.getChrom() == getChrom();
	}
	
	bool Interval::overlaps(const Interval& i) const {
		return i.getStart() < getEnd()
			and i.getEnd() > getStart()
			and i.getChrom() == getChrom();
	}

	Distance Interval::getLength() const {
		return getEnd() - getStart();
	}
	
	Coord Interval::getStartCoord() const {
		return Coord(getChrom(), getStart());
	}
	
	Coord Interval::getEndCoord() const {
		return Coord(getChrom(), getEnd());
	}

	std::ostream& operator<<(std::ostream& strm, const Interval& i) {
		return strm << i.getChrom()
					<< ':'
					<< i.getStart()
					<< '-'
					<< i.getEnd()
					<< i.getStrand();
	}

	bool operator<(const Interval& i1, const Interval& i2) {
		int chromComp = i1.getChrom().compare(i2.getChrom());
		return chromComp < 0
			or (chromComp == 0
				and (i1.getStart() < i2.getStart()
					 or (i1.getStart() == i2.getStart()
						 and i1.getEnd() < i2.getEnd())));
	}

	bool operator==(const Interval& i1, const Interval& i2) {
		return i1.getStart() == i2.getStart()
			and i1.getEnd() == i2.getEnd()
			and i1.getStrand() == i2.getStrand()
			and i1.getChrom() == i2.getChrom();
	}
	
	bool operator!=(const Interval& i1, const Interval& i2) {
		return not (i1 == i2);
	}

	bool operator<(const Coord& c, const Interval& i) {
		int chromComp = c.chrom.compare(i.getChrom());
		return chromComp < 0 or (chromComp == 0 and c.pos < i.getStart());
	}
	
} }
