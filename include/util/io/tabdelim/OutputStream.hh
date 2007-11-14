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

#ifndef __UTIL_IO_TABDELIM_OUTPUTSTREAM_HH__
#define __UTIL_IO_TABDELIM_OUTPUTSTREAM_HH__

#include <iosfwd>
#include <iterator>
#include <algorithm>

namespace util { namespace io { namespace tabdelim {

    class OutputStream {
	public:
		OutputStream(std::ostream& stream);

		template<typename C>
		OutputStream& operator<<(const C& fields);

		operator bool() const;
		bool operator!() const;

	private:
		std::ostream& stream;
	};

	inline OutputStream::OutputStream(std::ostream& stream) : stream(stream) {}
	inline OutputStream::operator bool() const { return stream; }
	inline bool OutputStream::operator!() const { return !stream; }

	template<typename C>
	OutputStream& OutputStream::operator<<(const C& fields) {
		std::copy(fields.begin(), fields.end(),
				  std::ostream_iterator<std::string>(stream, "\t"));
		stream << '\n';
		return *this;
	}

} } }

#endif // __UTIL_IO_TABDELIM_OUTPUTSTREAM_HH__
