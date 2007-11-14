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

#include <istream>
#include <sstream>
#include <stdexcept>

#include "bio/formats/repeatmasker/InputStream.hh"

namespace bio { namespace formats { namespace repeatmasker {

	InputStream::InputStream(std::istream& strm) : strm(strm) {
		skipLines(NUM_HEADER_LINES);
	}
	
	InputStream::operator bool() const { return strm; }
	
	bool InputStream::operator!() const { return !strm; }

	void InputStream::skipLines(size_t n) {
		for (size_t i = 0; i < n; ++i) {
			strm >> line;
		}
	}
	
	genome::Distance InputStream::parseLeft(const std::string& s) {
		if (s.size() < 3 or
			s[0] != LEFT_LEFT_DELIM or
			s[s.size() - 1] != LEFT_RIGHT_DELIM) {
			throw std::runtime_error("Invalid 'left' value: " + s);
		}
		std::istringstream leftStream(s.substr(1, s.size() - 2));
		genome::Distance left;
		leftStream >> left;
		return left;
	}
	
	InputStream& InputStream::operator>>(Record& rec) {
		strm >> line;
		if (not strm) { return *this; }
		
		std::istringstream lineStream(line);

		lineStream >> rec.score
				   >> rec.pctDivergence >> rec.pctDeleted >> rec.pctInserted
				   >> rec.queryName >> rec.queryStart >> rec.queryEnd;

		std::string queryLeft;
		std::string queryStrand;
		lineStream >> queryLeft >> queryStrand;
		rec.queryLeft = parseLeft(queryLeft);
		rec.queryStrand = (queryStrand == FORWARD ? '+' : '-');

		lineStream >> rec.repeatName >> rec.repeatClass;

		std::string repeatLeft;
		if (rec.queryStrand.isForward()) {
			lineStream >> rec.repeatStart >> rec.repeatEnd >> repeatLeft;
		} else {
			lineStream >> repeatLeft >> rec.repeatEnd >> rec.repeatStart;
		}
		rec.repeatLeft = parseLeft(repeatLeft);
		
		lineStream >> rec.id;

		if (not lineStream) {
			throw std::runtime_error("Invalid line:\n" + line);
		}

		std::string inHigherScoring;
		lineStream >> inHigherScoring;
		rec.isIncludedInHigherScoringMatch =
			(inHigherScoring == IN_HIGHER_SCORING);

		return *this;
	}

} } }
