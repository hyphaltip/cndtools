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

#include "bio/formats/agp/InputStream.hh"

namespace bio { namespace formats { namespace agp {

	void InputStream::stripComment(std::string& line) const {
		std::string::size_type pos = line.find('#');
		if (pos == std::string::npos) {
			return;
		} else {
			while (pos > 0 and std::isspace(line[pos - 1])) {
				--pos;
			}
			line.erase(pos);
		}
	}
	
	InputStream& InputStream::operator>>(Record& rec) {
		line.clear();
		while (line.empty() && strm) {
			strm >> line;
			stripComment(line);
		}
	
		if (not line.empty()) {
			line >> rec;
		}
	
		return *this;
	}

	InputStream& InputStream::operator>>(Record2& rec) {
		line.clear();
		while (line.empty() && strm) {
			strm >> line;
			stripComment(line);
		}
	
		if (line.empty()) { return *this; }
		
		std::istringstream lineStream(line);

		// Read first five fields (same for all record types)
		std::string type;
		lineStream >> rec.object
				   >> rec.objectStart
				   >> rec.objectEnd
				   >> rec.partNumber
				   >> type;

		rec.componentType = type.at(0);

		// Coordinates are 1-based closed
		--rec.objectStart;
		
		// Assumed to be on the forward strand
		rec.objectStrand = '+';

		// Read following fields according to whether this is a
		// gap record or not
		if (rec.isGap()) {
			std::string linkage;
			lineStream >> rec.gapLength >> rec.gapType >> linkage;
			rec.isLinkage = (linkage == LINKAGE_YES);
		} else {
			lineStream >> rec.componentID
					   >> rec.componentStart
					   >> rec.componentEnd
					   >> rec.componentStrand;

			// Coordinates are 1-base closed
			--rec.componentStart;
		}

		// Check to see parsing was successful
		if (not lineStream) {
			throw std::runtime_error("Invalid AGP line:\n" + line);
		}
	
		return *this;
	}

} } }
