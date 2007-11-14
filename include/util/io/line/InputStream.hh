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

#ifndef __UTIL_IO_LINE_INPUTSTREAM_HH__
#define __UTIL_IO_LINE_INPUTSTREAM_HH__

#include <iosfwd>
#include <vector>

#define IO_BUFFER_SIZE (4 * 1024)

namespace util { namespace io { namespace line {

	class InputStream {
	public:
		InputStream(std::istream& strm,
					const size_t bufferSize=IO_BUFFER_SIZE);
		InputStream& operator>>(std::string& line);
		operator bool() const;
		bool operator!() const;
				
	private:
		void fillBuffer();

		std::istream& strm;
		std::vector<char> buffer;
		std::vector<char>::iterator pos;
		bool atEnd;
	};
			
	//typedef util::stl::InputStreamIterator<InputStream> Iterator;

} } }

#endif // __UTIL_IO_LINE_INPUTSTREAM_HH__
