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

#ifndef __BIO_GENOME_INTERVALMAP_HH__
#define __BIO_GENOME_INTERVALMAP_HH__

#include <string>

#include "bio/genome/BasicInterval.hh"

namespace bio { namespace genome {

    class IntervalMap {
	public:
		virtual ~IntervalMap() {}

		virtual std::string getSourceChrom() const = 0;
		virtual Position getSourceStart() const = 0;
		virtual Position getSourceEnd() const = 0;
		virtual Strand getSourceStrand() const = 0;

		virtual std::string getTargetChrom() const = 0;
		virtual Position getTargetStart() const = 0;
		virtual Position getTargetEnd() const = 0;
		virtual Strand getTargetStrand() const = 0;

		virtual BasicInterval getSourceInterval() const;
		virtual BasicInterval getTargetInterval() const;
    };

} }

#endif // __BIO_GENOME_INTERVALMAP_HH__
