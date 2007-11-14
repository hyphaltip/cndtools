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

#ifndef __BIO_FORMATS_UCSC_CHAIN_RECORD_HH__
#define __BIO_FORMATS_UCSC_CHAIN_RECORD_HH__

#include <string>
#include <vector>

#include "bio/genome/Strand.hh"

namespace bio { namespace formats { namespace ucsc { namespace chain {

    class Record {
	public:
		struct Alignment {
			size_t size;
			size_t dt;
			size_t dq;

			void flip();
		};

		typedef std::vector<Alignment> AlignmentList;

		// Swap target and query
		void flip();
		
		long long score;
		std::string tName;
		size_t tSize;
		genome::Strand tStrand;
		size_t tStart;
		size_t tEnd;
		std::string qName;
		size_t qSize;
		genome::Strand qStrand;
		size_t qStart;
		size_t qEnd;
		std::string id;

		AlignmentList alignments;
	};

	inline void Record::Alignment::flip() { std::swap(dt, dq); }

} } } }

#endif // __BIO_FORMATS_UCSC_CHAIN_RECORD_HH__
