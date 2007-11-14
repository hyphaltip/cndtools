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

#include <vector>
#include <stdexcept>

#include "util/io.hh"

namespace util { namespace io {

	std::string readStream(std::istream& strm) {
		std::vector<char> buffer(IO_BUFFER_SIZE);
		std::string str;
		while (strm) {
			// Fill buffer
			strm.read(&buffer[0], buffer.size());

			// Append buffer to string
			str.append(buffer.begin(), buffer.begin() + strm.gcount());
		}
		return str;
	}

	void copy_stream(std::istream& in, std::ostream& out) {
		std::vector<char> buffer(IO_BUFFER_SIZE);
		while (in) {
			in.read(&buffer[0], buffer.size());
			out.write(&buffer[0], in.gcount());
		}
		// TODO: Check state of in stream and throw error if not reached EOF
	}
	
	std::string readStream(FILE* strm) {
		std::vector<char> buffer(IO_BUFFER_SIZE);
		std::string str;
		while (true) {
			// Fill buffer
			size_t count = fread(&buffer[0], sizeof(char), buffer.size(),
								 strm);

			// Check if we should stop
			if (count == 0) {
				break;
			}
				
			// Append buffer to string
			str.append(buffer.begin(), buffer.begin() + count);
		}
		return str;
	}

} }
