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

#include "bio/formats/agp/Record2.hh"

namespace bio { namespace formats { namespace agp {

	bool Record2::isGap() const {
		return componentType == GAP;
	}

	// Interval interface
	std::string Record2::getChrom() const { return object; }
	genome::Position Record2::getStart() const { return objectStart; }
	genome::Position Record2::getEnd() const { return objectEnd; }
	genome::Strand Record2::getStrand() const { return objectStrand; }
	void Record2::setChrom(const std::string& chrom) {
		object = chrom;
	}
	void Record2::setStart(const genome::Position start) {
		objectStart = start;
	}
	void Record2::setEnd(const genome::Position end) {
		objectEnd = end;
	}
	void Record2::setStrand(const genome::Strand strand) {
		objectStrand = strand;
	}
	
	// Interval map interface
	std::string Record2::getSourceChrom() const { return componentID; }
	genome::Position Record2::getSourceStart() const { return componentStart; }
	genome::Position Record2::getSourceEnd() const { return componentEnd; }
	genome::Strand Record2::getSourceStrand() const { return componentStrand; }
	
	std::string Record2::getTargetChrom() const { return object; }
	genome::Position Record2::getTargetStart() const { return objectStart; }
	genome::Position Record2::getTargetEnd() const { return objectEnd; }
	genome::Strand Record2::getTargetStrand() const { return objectStrand; }

	void Record2::setSourceChrom(const std::string& chrom) {
		componentID = chrom;
	}
	
	void Record2::setSourceStart(const genome::Position start) {
		componentStart = start;
	}
	
	void Record2::setSourceEnd(const genome::Position end) {
		componentEnd = end;
	}
	
	void Record2::setSourceStrand(const genome::Strand strand) {
		componentStrand = strand;
	}
	
	void Record2::setTargetChrom(const std::string& chrom) {
		object = chrom;
	}
	
	void Record2::setTargetStart(const genome::Position start) {
		objectStart = start;
	}
	
	void Record2::setTargetEnd(const genome::Position end) {
		objectEnd = end;
	}
	
	void Record2::setTargetStrand(const genome::Strand strand) {
		objectStrand = strand;
	}
	
} } }
