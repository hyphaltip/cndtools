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

#include <cassert>

#include "bio/formats/clustal/Record.hh"

namespace bio { namespace formats { namespace clustal {

	size_t Record::getNumCols() const {
		assert(not sequences.empty());
		return sequences.front().size();
	}
	
	size_t Record::getNumSeqs() const {
		return sequences.size();
	}
	
	std::string Record::getSeq(size_t seqNum) const {
		return sequences[seqNum];
	}
	   
	std::string Record::getSubstring(size_t seqNum,
									 size_t colStart,
									 size_t len) const {
		return getSeq(seqNum).substr(colStart, len);
	}
	
	char Record::getChar(size_t seqNum, size_t colNum) const {
		return getSeq(seqNum)[colNum];
	}

	std::string Record::getName(size_t seqNum) const {
		return names[seqNum];
	}

} } }
