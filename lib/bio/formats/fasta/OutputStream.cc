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

#include <ostream>

#include "bio/formats/fasta/OutputStream.hh"

namespace bio { namespace formats { namespace fasta {

	OutputStream::OutputStream(std::ostream& strm)
		: strm(strm),
		  lineWidth(DEFAULT_LINE_WIDTH) {
	}
		
	void OutputStream::setLineWidth(size_t lineWidth) {
		this->lineWidth = lineWidth;
	}
	
	OutputStream& OutputStream::operator<<(const Record& r) {
		// Output title line
		strm << TITLE_LINE_PREFIX << r.title << '\n';

		// Output sequence wrapped at desired width
		for (size_t start = 0; start < r.sequence.size(); start += lineWidth) {
			strm << r.sequence.substr(start, lineWidth) << '\n';
		}

		return *this;
	}

	OutputStream& OutputStream::operator<<(const alignment::NamedMultipleAlignment& a) {
		for (size_t i = 0; i < a.getNumSeqs(); ++i) {
			(*this) << Record(a.getName(i), a.getSeq(i));
		}
		return *this;
	}
		
} } }
