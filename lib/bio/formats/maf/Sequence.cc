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

#include "bio/formats/maf/Sequence.hh"
#include <algorithm>

namespace bio { namespace formats { namespace maf {

	std::string Sequence::getChrom() const {
		return src;
	}
	
	genome::Position Sequence::getStart() const {
		if (strand.isForward()) {
			return start;
		} else {
			return srcSize - start - size;
		}
	}
			   
	genome::Position Sequence::getEnd() const {
		if (strand.isForward()) {
			return start + size;
		} else {
			return srcSize - start;
		}
	}
		
	genome::Strand Sequence::getStrand() const {
		return strand;
	}

	void Sequence::setSequence(const std::string& s) {
		text = s;
		size = s.size() - std::count(s.begin(), s.end(), '-');
	}
	
	void Sequence::setChromSize(const genome::Distance size) {
		genome::Position absStart = getStart();
		srcSize = size;
		setStart(absStart);
	}
	
	void Sequence::setChrom(const std::string& chrom) {
		src = chrom;
	}

	void Sequence::setStart(const genome::Position start) {
		if (strand.isForward()) {
			this->start = start;
		} else {
			this->start = srcSize - start - size;
		}
	}

	void Sequence::setEnd(const genome::Position end) {
		genome::Position absStart = end - size;
		setStart(absStart);
	}

	void Sequence::setStrand(const genome::Strand strand) {
		genome::Position absStart = getStart();
		this->strand = strand;
		setStart(absStart);
	}

	std::pair<std::string, std::string> Sequence::getGenomeAndChrom() const {
		std::string::size_type dot_pos = src.find('.');
		return std::make_pair(src.substr(0, dot_pos), src.substr(dot_pos + 1));
	}
	
} } }
