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

#include <sstream>

#include "bio/formats/ucsc/chain/InputStream.hh"
#include "util/string.hh"

namespace bio { namespace formats { namespace ucsc { namespace chain {

	using util::string::startsWith;
	
	void InputStream::readToNextNonCommentLine() {
		const std::string COMMENT_PREFIX = "#";
		do {
			stream >> line;
		} while (stream and startsWith(line, COMMENT_PREFIX));
	}
	
	InputStream& InputStream::operator>>(Record& r) {
		const std::string CHAIN_HEADER_PREFIX = "chain ";

		readToNextNonCommentLine();
		if (not stream) { return *this; }
		
		if (not startsWith(line, CHAIN_HEADER_PREFIX)) {
			throw std::runtime_error("Invalid chain header line:\n" + line);
		}

		std::string header = line.substr(CHAIN_HEADER_PREFIX.size());

		std::istringstream headerStream(header);
		headerStream >> r.score
					 >> r.tName >> r.tSize >> r.tStrand >> r.tStart >> r.tEnd
					 >> r.qName >> r.qSize >> r.qStrand >> r.qStart >> r.qEnd
					 >> r.id;
		
		if (not headerStream) {
			throw std::runtime_error("Invalid chain header line:\n" + line);
		}

		r.alignments.clear();
		
		Record::Alignment a;

		while (true) {
			readToNextNonCommentLine();
			if (not stream) {
				throw std::runtime_error("EOF encountered before chain end");
			}
			if (line.empty()) { break; }
			std::istringstream alignStream(line);
			alignStream >> a.size >> a.dt >> a.dq;
			r.alignments.push_back(a);
		}

		// Set last alignment block's gaps to zero
		if (r.alignments.size() > 0) {
			r.alignments.back().dt = 0;
			r.alignments.back().dq = 0;
		}

		return *this;
	}
	
} } } }
