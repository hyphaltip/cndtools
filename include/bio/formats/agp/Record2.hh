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

#ifndef __BIO_FORMATS_AGP_RECORD2_HH__
#define __BIO_FORMATS_AGP_RECORD2_HH__

#include <string>
#include <iosfwd>
#include <vector>

#include "bio/genome/MutableInterval.hh"
#include "bio/genome/MutableIntervalMap.hh"
#include "bio/formats/agp/Constants.hh"

// Facilities for processing AGP formatted records
// See http://genome.ucsc.edu/goldenPath/datorg.html for information

namespace bio { namespace formats { namespace agp {

	class Record2 : public genome::MutableInterval,
					public genome::MutableIntervalMap,
					public Constants {
	public:
		std::string object;
		genome::Position objectStart;
		genome::Position objectEnd;
		genome::Strand objectStrand;

		size_t partNumber;
		char componentType;

		std::string componentID;
		genome::Position componentStart;
		genome::Position componentEnd;
		genome::Strand componentStrand;

		genome::Distance gapLength;
		std::string gapType;
		bool isLinkage;

		bool isGap() const;
		
		// Interval interface
		std::string getChrom() const;
		genome::Position getStart() const;
		genome::Position getEnd() const;
		genome::Strand getStrand() const;
		void setChrom(const std::string& chrom);
		void setStart(const genome::Position start);
		void setEnd(const genome::Position end);
		void setStrand(const genome::Strand strand);

		// Interval map interface
		std::string getSourceChrom() const;
		genome::Position getSourceStart() const;
		genome::Position getSourceEnd() const;
		genome::Strand getSourceStrand() const;

		std::string getTargetChrom() const;
		genome::Position getTargetStart() const;
		genome::Position getTargetEnd() const;
		genome::Strand getTargetStrand() const;

		void setSourceChrom(const std::string& chrom);
		void setSourceStart(const genome::Position start);
		void setSourceEnd(const genome::Position end);
		void setSourceStrand(const genome::Strand strand);

		void setTargetChrom(const std::string& chrom);
		void setTargetStart(const genome::Position start);
		void setTargetEnd(const genome::Position end);
		void setTargetStrand(const genome::Strand strand);
	};


} } }

#endif // __BIO_FORMATS_AGP_RECORD2_HH__
