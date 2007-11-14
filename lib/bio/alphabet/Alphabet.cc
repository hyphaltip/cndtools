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

#include "bio/alphabet/Alphabet.hh"

namespace bio { namespace alphabet {
	
	bool Alphabet::isOn(const std::string& s) const {
		for (std::string::const_iterator it = s.begin(); it != s.end(); ++it) {
			if (!isMember(*it)) {
				return false;
			}
		}
		return true;
	}

	std::string Alphabet::encode(const std::string& s) const {
		std::string e(s);
		for (std::string::iterator i = e.begin(); i != e.end(); ++i) {
			*i = encode(*i);
		}
		return e;
	}

	std::string Alphabet::decode(const std::string& s) const {
		std::string d(s);
		for (std::string::iterator i = d.begin(); i != d.end(); ++i) {
			*i = decode(*i);
		}
		return d;
	}
		
} }
