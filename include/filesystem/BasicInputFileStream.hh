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

#ifndef __FILESYSTEM_BASICINPUTFILESTREAM_HH__
#define __FILESYSTEM_BASICINPUTFILESTREAM_HH__

#include <fstream>
#include <stdexcept>

#include "filesystem/Path.hh"

namespace filesystem {

    template <class charT, class traits = std::char_traits<charT> >
    class BasicInputFileStream : public std::basic_ifstream<charT,traits> {
    public:
		BasicInputFileStream() {}
		explicit BasicInputFileStream(const Path& path,
									  std::ios_base::openmode mode = std::ios_base::in)
			: std::basic_ifstream<charT, traits>(path.toString().c_str(),
												 mode) {
			if (not *this) {
				throw std::runtime_error("Could not open input file: "
										 + path.toString());
			}
		}

		void open(const Path& path,
				  std::ios_base::openmode mode = std::ios_base::in) {
			std::basic_ifstream<charT, traits>::open(path.toString().c_str(),
													 mode);
			if (not *this) {
				throw std::runtime_error("Could not open input file: "
										 + path.toString());
			}
		}

		virtual ~BasicInputFileStream() {}
    };

}

#endif // __FILESYSTEM_BASICINPUTFILESTREAM_HH__
