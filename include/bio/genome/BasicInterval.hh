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

#ifndef __BIO_GENOME_BASICINTERVAL_HH__
#define __BIO_GENOME_BASICINTERVAL_HH__

#include "bio/genome/MutableInterval.hh"

namespace bio { namespace genome {

    class BasicInterval : public MutableInterval {
	public:
		static BasicInterval between(const Interval& i1, const Interval& i2);
		
		virtual std::string getChrom() const;
		virtual Position getStart() const;
		virtual Position getEnd() const;
		virtual Strand getStrand() const;

		virtual void setChrom(const std::string& chrom);
		virtual void setStart(const Position start);
		virtual void setEnd(const Position end);
		virtual void setStrand(const Strand strand);

		BasicInterval(const Interval& i);
		
		BasicInterval(const std::string& chrom = "",
					  const Position start = 0,
					  const Position end = 0,
					  const Strand strand = '+');

		BasicInterval(const Coord& start,
					  const Coord& end,
					  const Strand strand = '+');

		BasicInterval(const Coord& c,
					  const Strand strand = '+');
		
		BasicInterval operator&(const Interval& i) const;
		BasicInterval operator|(const Interval& i) const;
		BasicInterval& operator&=(const Interval& i);
		BasicInterval& operator|=(const Interval& i);
		
	private:
		std::string chrom;
		Position start;
		Position end;
		Strand strand;
	};
	
} }

#endif // __BIO_GENOME_BASICINTERVAL_HH__
