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

#ifndef __BIO_GENOME_INTERVAL_HH__
#define __BIO_GENOME_INTERVAL_HH__

#include <string>
#include <iosfwd>

#include "bio/genome/Coord.hh"
#include "bio/genome/Strand.hh"

namespace bio { namespace genome {

	class Interval {
	public:
		virtual ~Interval();
		
		virtual std::string getChrom() const = 0;
		virtual Position getStart() const = 0;
		virtual Position getEnd() const = 0;
		virtual Strand getStrand() const = 0;

		virtual Coord getStartCoord() const;
		virtual Coord getEndCoord() const;
		virtual Distance getLength() const;
		
		virtual bool contains(const Coord& c) const;
		virtual bool contains(const Interval& i) const;
		virtual bool overlaps(const Interval& i) const;
	};
	
	std::ostream& operator<<(std::ostream& strm, const Interval& i);
	
	bool operator<(const Interval& i1, const Interval& i2);
	bool operator==(const Interval& i1, const Interval& i2);
	bool operator!=(const Interval& i1, const Interval& i2);
	
	bool operator<(const Coord& c, const Interval& i);
	
} }

#endif // __BIO_GENOME_INTERVAL_HH__
