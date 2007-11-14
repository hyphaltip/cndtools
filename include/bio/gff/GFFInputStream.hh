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

#ifndef __GFF_INPUT_STREAM_HH__
#define __GFF_INPUT_STREAM_HH__

#include <string>
#include <iosfwd>

#include "util/io/line/InputStream.hh"
#include "bio/gff/GFFRecord.hh"

namespace bio { namespace gff {

	class GFFInputStream {
	public:
		typedef GFFRecord ValueType;
	private:
		util::io::line::InputStream strm;
		std::string line;

	public:
		GFFInputStream(std::istream& strm,
					   const size_t bufferSize=IO_BUFFER_SIZE);

		GFFInputStream& operator>>(GFFRecord& rec);
		operator bool();
		bool operator!();
	};

	inline GFFInputStream::GFFInputStream(std::istream& strm,
										  const size_t bufferSize) :
		strm(strm, bufferSize) {
	}
	
	inline GFFInputStream::operator bool() {
		return strm;
	}
	
	inline bool GFFInputStream::operator!() {
		return !strm;
	}

} }

#endif // __GFF_INPUT_STREAM_HH__
