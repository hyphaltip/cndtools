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

#ifndef __BIO_GENOME_MUTABLEINTERVAL_HH__
#define __BIO_GENOME_MUTABLEINTERVAL_HH__

#include "bio/genome/Interval.hh"

namespace bio { namespace genome {

    class MutableInterval : public Interval {
	public:
		virtual void setChrom(const std::string& chrom) = 0;
		virtual void setStart(const Position start) = 0;
		virtual void setEnd(const Position end) = 0;
		virtual void setStrand(const Strand strand) = 0;

		virtual void setStartCoord(const Coord& c);
		virtual void setEndCoord(const Coord& c);
		virtual void setInterval(const Interval& i);

		// Change this interval to the opposite strand
		virtual void flip();
		
		virtual MutableInterval& operator=(const Interval& i);
    };

	std::istream& operator>>(std::istream& strm, MutableInterval& i);

} }

#endif // __BIO_GENOME_MUTABLEINTERVAL_HH__
