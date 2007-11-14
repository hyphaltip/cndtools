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

#ifndef __BIO_FORMATS_REPEATMASKER_INPUTSTREAM_HH__
#define __BIO_FORMATS_REPEATMASKER_INPUTSTREAM_HH__

#include <string>
#include <iosfwd>

#include "util/io/line/InputStream.hh"
#include "bio/formats/repeatmasker/Constants.hh"
#include "bio/formats/repeatmasker/Record.hh"

namespace bio { namespace formats { namespace repeatmasker {

    class InputStream : private Constants {
	public:
		InputStream(std::istream& strm);

		InputStream& operator>>(Record& rec);
		operator bool() const;
		bool operator!() const;
		
	private:
		static genome::Distance parseLeft(const std::string& s);
		void skipLines(size_t n);

		util::io::line::InputStream strm;
		std::string line;
	};
	
} } }

#endif // __BIO_FORMATS_REPEATMASKER_INPUTSTREAM_HH__
