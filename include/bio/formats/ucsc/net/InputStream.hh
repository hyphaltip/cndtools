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

#ifndef __BIO_FORMATS_UCSC_NET_INPUTSTREAM_HH__
#define __BIO_FORMATS_UCSC_NET_INPUTSTREAM_HH__

#include "util/io/line/InputStream.hh"
#include "bio/formats/ucsc/net/Record.hh"

namespace bio { namespace formats { namespace ucsc { namespace net {

	class InputStream {
	public:
		InputStream(std::istream& stream,
					const size_t bufferSize=IO_BUFFER_SIZE);

		InputStream& operator>>(Record& rec);
		operator bool();
		bool operator!();
		
	private:
		util::io::line::InputStream stream;
		std::string line;
		std::string tChrom;
		Record::Distance tChromSize;
	};

	inline InputStream::InputStream(std::istream& stream,
									const size_t bufferSize) :
		stream(stream, bufferSize) {
	}
	
	inline InputStream::operator bool() {
		return stream;
	}
	
	inline bool InputStream::operator!() {
		return !stream;
	}

} } } }

#endif // __BIO_FORMATS_UCSC_NET_INPUTSTREAM_HH__
