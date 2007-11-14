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

#ifndef __BIO_FORMATS_UCSC_NET_RECORD_HH__
#define __BIO_FORMATS_UCSC_NET_RECORD_HH__

#include "bio/genome/Interval.hh"
#include "util/stl.hh"

namespace bio { namespace formats { namespace ucsc { namespace net {

	class Record {
	public:
		typedef genome::Distance Distance;
		typedef genome::Strand Orientation;
		typedef util::stl::hash_map<std::string, std::string> AttributeMap;

		enum RecClass {
			FILL, // A portion of a chain.
			GAP,
		};

		// Fixed fields
		RecClass recClass;   // Either fill or gap. Fill refers to a
						     // portion of a chain.
		std::string tChrom;	 // Chromosome name (target species)
		Distance tChromSize; // Size of chromosome (target species)
		Distance tStart;     // Start in chromosome (target species)
		Distance tSize;      // Size (target species)
		std::string qChrom;  // Chromsome name (query species)
		Distance qStart;     // Start in chromsome (query species)
		Distance qSize;      // Size (query species)
		Orientation orientation;  // Relative orientation between target and
	                              // query species.
		unsigned int level;

		// * id -- ID of associated chain (gapped alignment), if any.
		// * score -- Score of associated chain.
		// * ali -- Number of bases in alignments in chain.
		// * qFar -- For fill that is on the same chromosome as parent, how far fill is from position predicted by parent. This helps determine if a rearrangement is local or if a duplication is tandem.
		// * qOver -- Number of bases overlapping with parent gap on query side. Generally, this will be near zero, except for inverts.
		// * qDup -- Number of bases in query region that are used twice or more in net. This helps distinguish between a rearrangement and a duplication.
		// * type -- One of the following values:
		//    o top -- Chain is top-level, not a gap filler
		//    o syn -- Chain is on same chromosome and in same direction as parent
		//    o inv -- Chain is on same chromosome on opposite direction from parent
		//    o nonSyn -- Chain is on a different chromosome from parent 
		// * tN -- Number of unsequenced bases (Ns) on target side
		// * qN -- Number of unsequenced bases on query side
		// * tR -- Number of bases in RepeatMasker masked repeats on target.
		// * qR -- Number of bases in RepeatMasker masked repeats on query.
		// * tNewR -- Bases in lineage-specific repeats on target.
		// * qNewR -- Bases in lineage-specific repeats on query.
		// * tOldR -- Bases in repeats predating split on target.
		// * qOldR -- Bases in repeats predating split on query.
		// * tTrf -- Bases in trf (Tandem Repeat Finder) repeats on target.
		// * qTrf -- Bases in trf repeats on query. 		
		AttributeMap attributes;
    };

} } } }

#endif // __BIO_FORMATS_UCSC_NET_RECORD_HH__
