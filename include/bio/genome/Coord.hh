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

#ifndef __BIO_GENOME_COORD_HH__
#define __BIO_GENOME_COORD_HH__

#include <string>

namespace bio { namespace genome {

	typedef long long Position;
	typedef long long Distance;

	class Coord {
	public:
		std::string chrom;
		Position pos;

		Coord(const std::string& chrom="", const Position pos=0);

		Coord& operator+=(Distance d);
		Coord operator+(Distance d) const;
		Coord& operator-=(Distance d);
		Coord operator-(Distance d) const;
	};
	
	bool operator<(const Coord& c1, const Coord& c2);
	bool operator<=(const Coord& c1, const Coord& c2);
	bool operator==(const Coord& c1, const Coord& c2);
	Distance operator-(const Coord& c1, const Coord& c2);
	
	inline Coord::Coord(const std::string& chrom, const Position pos)
		: chrom(chrom), pos(pos) {
	}

	inline bool operator<(const Coord& c1, const Coord& c2) {
		int chromComp = c1.chrom.compare(c2.chrom);
		return chromComp < 0 or (chromComp == 0 and c1.pos < c2.pos);
	}

	inline bool operator<=(const Coord& c1, const Coord& c2) {
		int chromComp = c1.chrom.compare(c2.chrom);
		return chromComp < 0 or (chromComp == 0 and c1.pos <= c2.pos);
	}

	inline bool operator==(const Coord& c1, const Coord& c2) {
		return c1.pos == c2.pos and c1.chrom == c2.chrom;
	}

	inline Distance operator-(const Coord& c1, const Coord& c2) {
		return c1.pos - c2.pos;
	}

	inline Coord& Coord::operator+=(Distance d) {
		pos += d;
		return *this;
	}
	
	inline Coord Coord::operator+(Distance d) const {
		return Coord(*this) += d;
	}

	inline Coord& Coord::operator-=(Distance d) {
		pos -= d;
		return *this;
	}
	
	inline Coord Coord::operator-(Distance d) const {
		return Coord(*this) -= d;
	}
	
} }

#endif // __BIO_GENOME_COORD_HH__
