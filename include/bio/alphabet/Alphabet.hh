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

#ifndef __ALPHABET_HH__
#define __ALPHABET_HH__

#include <string>

namespace bio {

namespace alphabet {

	class Alphabet {
	public:
		virtual ~Alphabet() {}
		virtual unsigned char getSize() const = 0;
		virtual bool isMember(const char c) const = 0;
		bool isOn(const std::string& s) const;
		virtual unsigned char encode(const char c) const = 0;
		virtual std::string encode(const std::string& s) const;
		virtual char decode(const unsigned char c) const = 0;
		virtual std::string decode(const std::string& s) const;
	};

};

};

#endif // __ALPHABET_HH__
