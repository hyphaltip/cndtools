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

#ifndef __BIO_FORMATS_AGP_INPUTSTREAM_HH__
#define __BIO_FORMATS_AGP_INPUTSTREAM_HH__

#include <string>
#include <iosfwd>

#include "util/io/line/InputStream.hh"
#include "util/string.hh"
#include "bio/formats/agp/Record.hh"
#include "bio/formats/agp/Record2.hh"
#include "bio/formats/agp/Constants.hh"

namespace bio { namespace formats { namespace agp {

	class InputStream : private Constants {
	public:
		InputStream(std::istream& strm,
					const size_t bufferSize=IO_BUFFER_SIZE);

		InputStream& operator>>(Record& rec);
		InputStream& operator>>(Record2& rec);
		operator bool();
		bool operator!();
		
	private:
		void stripComment(std::string& line) const;
		
		util::io::line::InputStream strm;
		std::string line;
	};

	inline InputStream::InputStream(std::istream& strm,
									const size_t bufferSize) :
		strm(strm, bufferSize) {
	}
	
	inline InputStream::operator bool() {
		return strm;
	}
	
	inline bool InputStream::operator!() {
		return !strm;
	}

} } }

#endif // __BIO_FORMATS_AGP_INPUTSTREAM_HH__
