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

#include "bio/formats/repeatmasker/Constants.hh"

namespace bio { namespace formats { namespace repeatmasker {

	const std::string Constants::HEADER =
		"   SW  perc perc perc  query     position in query            matching       repeat         position in  repeat\n"
		"score  div. del. ins.  sequence    begin     end    (left)   repeat         class/family    begin  end (left)   ID\n"
		"\n";
	
	const size_t Constants::NUM_HEADER_LINES = 3;
	const char Constants::LEFT_LEFT_DELIM = '(';
	const char Constants::LEFT_RIGHT_DELIM = ')';
	const std::string Constants::FORWARD = "+";
	const std::string Constants::REVERSE = "C";
	const std::string Constants::IN_HIGHER_SCORING = "*";

} } }
