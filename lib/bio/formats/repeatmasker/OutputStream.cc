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

#include "bio/formats/repeatmasker/OutputStream.hh"

namespace bio { namespace formats { namespace repeatmasker {

	OutputStream::OutputStream(std::ostream& strm) : strm(strm) {
		strm << HEADER;
	}
	
	OutputStream& OutputStream::operator<<(const Record& r) {
		strm << r.score << '\t'
			 << r.pctDivergence << '\t'
			 << r.pctDeleted << '\t'
			 << r.pctInserted << '\t'
			 << r.queryName << '\t'
			 << r.queryStart << '\t'
			 << r.queryEnd << '\t'
			 << LEFT_LEFT_DELIM << r.queryLeft << LEFT_RIGHT_DELIM << '\t'
			 << (r.queryStrand.isForward() ? FORWARD : REVERSE) << '\t'
			 << r.repeatName << '\t'
			 << r.repeatClass << '\t';
		if (r.queryStrand.isForward()) {
			strm << r.repeatStart << '\t'
				 << r.repeatEnd << '\t'
				 << LEFT_LEFT_DELIM << r.queryLeft << LEFT_RIGHT_DELIM << '\t';
		} else {
			strm << LEFT_LEFT_DELIM << r.queryLeft << LEFT_RIGHT_DELIM << '\t'
				 << r.repeatEnd << '\t'
				 << r.repeatStart << '\t';
		}
		strm << r.id;
		if (r.isIncludedInHigherScoringMatch) {
			strm << '\t' << IN_HIGHER_SCORING;
		}
		strm << '\n';

		return *this;
	}
	
	OutputStream::operator bool() const { return strm; }
	bool OutputStream::operator!() const { return not strm; }
	
} } }
