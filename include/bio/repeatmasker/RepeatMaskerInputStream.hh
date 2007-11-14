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

#ifndef __BIO_REPEATMASKER_REPEATMASKERINPUTSTREAM_HH__
#define __BIO_REPEATMASKER_REPEATMASKERINPUTSTREAM_HH__

#include <string>
#include <iosfwd>

#include "util/io/line/InputStream.hh"
#include "bio/repeatmasker/RepeatMaskerRecord.hh"

namespace bio { namespace repeatmasker {

    class RepeatMaskerInputStream {
	public:
		RepeatMaskerInputStream(std::istream& strm,
								const size_t bufferSize=(IO_BUFFER_SIZE << 3));

		RepeatMaskerInputStream& operator>>(RepeatMaskerRecord& rec);
		operator bool() const;
		bool operator!() const;
		
		typedef RepeatMaskerRecord ValueType;
		
	private:
		util::io::line::InputStream strm;
		std::string line;
	};
	
	inline RepeatMaskerInputStream::operator bool() const {
		return strm;
	}
	
	inline bool RepeatMaskerInputStream::operator!() const {
		return !strm;
	}
		
} }

#endif // __BIO_REPEATMASKER_REPEATMASKERINPUTSTREAM_HH__
