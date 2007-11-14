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

#ifndef __BIO_GENOME_MUTABLEINTERVALMAP_HH__
#define __BIO_GENOME_MUTABLEINTERVALMAP_HH__

#include "bio/genome/IntervalMap.hh"

namespace bio { namespace genome {

    class MutableIntervalMap : public IntervalMap {
		virtual void setSourceChrom(const std::string& chrom) = 0;
		virtual void setSourceStart(const Position start) = 0;
		virtual void setSourceEnd(const Position end) = 0;
		virtual void setSourceStrand(const Strand strand) = 0;

		virtual void setTargetChrom(const std::string& chrom) = 0;
		virtual void setTargetStart(const Position start) = 0;
		virtual void setTargetEnd(const Position end) = 0;
		virtual void setTargetStrand(const Strand strand) = 0;

		virtual void setSourceInterval(const Interval& i);
		virtual void setTargetInterval(const Interval& i);
    };

} }

#endif // __BIO_GENOME_MUTABLEINTERVALMAP_HH__
