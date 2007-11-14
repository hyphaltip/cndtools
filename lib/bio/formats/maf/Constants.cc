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

#include "bio/formats/maf/Constants.hh"

namespace bio { namespace formats { namespace maf {

	const std::string Constants::VERSION = "1";
	const std::string Constants::ALIGNMENT_LINE_PREFIX = "a";
	const std::string Constants::SEQ_LINE_PREFIX = "s";
	const std::string Constants::HEADER_LINE_PREFIX = "##maf";
	const std::string Constants::METADETA_LINE_PREFIX = "##";
	const std::string Constants::COMMENT_LINE_PREFIX = "#";
	
} } }
