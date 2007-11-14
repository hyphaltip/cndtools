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

#include "bio/formats/fasta/InputStream.hh"
#include "util/string.hh"
using util::string::IsWhitespace;

namespace bio { namespace formats { namespace fasta {
	
	InputStream::operator bool() {
		return !atEnd;
	}

	bool InputStream::operator!() {
		return atEnd;
	}

	void InputStream::setKeepWhitespace(bool b) {
 		keepWhitespace = b;
	}

	InputStream::InputStream(std::istream& strm,
							 const size_t bufferSize)
		: strm(strm),
		  buffer(bufferSize),
		  pos(buffer.end()),
		  atEnd(false),
		  keepWhitespace(false)
	{
		// Fill buffers until first '>' character is found or EOF
		while (strm && pos == buffer.end()) {
			fillBuffer();
			pos = std::find(buffer.begin(), buffer.end(), TITLE_LINE_PREFIX);
		}

		// If no > was found, mark stream as processed
		if (pos == buffer.end()) {
			atEnd = true;
		}
	}

	void InputStream::fillBuffer() {
		// Fill buffer
		strm.read(&buffer[0], buffer.size());

		// If buffer is not completely filled, resize to fit
		if (static_cast<size_t>(strm.gcount()) < buffer.size()) {
			buffer.resize(strm.gcount());
		}
	}
	
	InputStream& InputStream::operator>>(Record& rec) {
		// Return immediately if the stream has been processed
		if (pos == buffer.end()) {
			atEnd = true;
			return *this;
		}

		rec.title.clear();
		
		// Skip over '>' character
		std::vector<char>::iterator start = pos + 1;

		// Find newline char
		pos = std::find(start, buffer.end(), '\n');

		// Fill buffers until end of title is found
		while (strm && pos == buffer.end()) {
			rec.title.append(start, pos);
			fillBuffer();
			start = buffer.begin();
			pos = std::find(start, buffer.end(), '\n');
		}
		rec.title.append(start, pos);

		// Remove trailing whitespace from the title
		rec.title = util::string::stripRight(rec.title);

		rec.sequence.clear();

		// Skip over newline
		start = pos + 1;

		// Find the next title line
		pos = std::find(start, buffer.end(), TITLE_LINE_PREFIX);
		
		// Fill buffers until end of sequence is found
		while (strm && pos == buffer.end()) {
			if (keepWhitespace) {
				rec.sequence.append(start, pos);
			} else {
				rec.sequence.append(start,
									std::remove_if(start, pos, IsWhitespace()));
			}
			fillBuffer();
			start = buffer.begin();
			pos = std::find(start, buffer.end(), TITLE_LINE_PREFIX);
		}

		if (keepWhitespace) {
			rec.sequence.append(start, pos);
		} else {
			rec.sequence.append(start,
								std::remove_if(start, pos, IsWhitespace()));
		}
		
		return *this;
	}

	InputStream& InputStream::operator>>(alignment::BasicNamedMultipleAlignment& a) {
		a.clear();
		Record r;
		while ((*this) >> r) {
			a.addSeq(r.sequence, r.title);
		}
		return *this;
	}

} } }
