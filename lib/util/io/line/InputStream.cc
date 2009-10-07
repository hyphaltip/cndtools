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

#include "util/io/line/InputStream.hh"

#include <istream>
#include <algorithm>

namespace util { namespace io { namespace line {

	InputStream::InputStream(std::istream& strm, const size_t bufferSize)
		: strm(strm),
		  buffer(bufferSize),
		  pos(buffer.end()),
		  atEnd(false) {
		// Fill buffer and initialize position to start
		fillBuffer();
		pos = buffer.begin();

		// Mark stream as processed if input stream was empty
		if (pos == buffer.end()) {
			atEnd = true;
		}
	}

	InputStream::operator bool() const {
		return !atEnd;
	}
			
	bool InputStream::operator!() const {
		return atEnd;
	}
	
	void InputStream::fillBuffer() {
		// Fill buffer
		strm.read(&buffer[0], buffer.size());

		// If buffer is not completely filled, resize to fit
		if (static_cast<size_t>(strm.gcount()) < buffer.size()) {
			buffer.resize(strm.gcount());
		}
	}
	
	InputStream& InputStream::operator>>(std::string& line) {
		// Return immediately if the stream has been processed
		if (pos == buffer.end()) {
			atEnd = true;
			return *this;
		}

		// Empty the line string
		line.clear();

		// Mark beginning of line
		std::vector<char>::iterator start = pos;
		
		// Find newline char
		pos = std::find(start, buffer.end(), '\n');

		// Fill buffers until end of line is found
		while (strm && pos == buffer.end()) {
			line.append(start, pos);
			fillBuffer();
			start = buffer.begin();
			pos = std::find(start, buffer.end(), '\n');
		}
		line.append(start, pos);

		// Unless we hit the end of the file, skip over the next newline
		if (pos != buffer.end()) {
			++pos;
			if (pos == buffer.end() && strm) {
				fillBuffer();
				pos = buffer.begin();
			}
		}
		
		return *this;
	}

} } }
