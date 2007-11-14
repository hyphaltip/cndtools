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

#include "bio/formats/agp/OutputStream.hh"

namespace bio { namespace formats { namespace agp {

	OutputStream::OutputStream(std::ostream& strm) : strm(strm) {
	}
		
	OutputStream::operator bool() const { return strm; }
	bool OutputStream::operator!() const { return not strm; }

	OutputStream& OutputStream::operator<<(const Record2& rec) {
		strm << rec.object << '\t'
			 << rec.objectStart + 1 << '\t'
			 << rec.objectEnd << '\t'
			 << rec.partNumber << '\t'
			 << rec.componentType << '\t';
		if (rec.isGap()) {
			strm << rec.gapLength << '\t'
				 << rec.gapType << '\t'
				 << (rec.isLinkage ? LINKAGE_YES : LINKAGE_NO) << '\n';
		} else {
			strm << rec.componentID << '\t'
				 << rec.componentStart + 1 << '\t'
				 << rec.componentEnd << '\t'
				 << (rec.objectStrand.isForward() ?
					 rec.componentStrand :
					 rec.componentStrand.opposite())
				 << '\n';
		}
	
		return *this;
	}

} } }
