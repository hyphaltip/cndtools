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

#include "bio/formats/agp/Constants.hh"

namespace bio { namespace formats { namespace agp {

	// Component types
	const char Constants::ACTIVE_FINISHING = 'A';
	const char Constants::DRAFT_HTG = 'D';
	const char Constants::FINISHED_HTG = 'F';
	const char Constants::WHOLE_GENOME_FINISHING = 'G';
	const char Constants::PRE_DRAFT = 'P';
	const char Constants::GAP = 'N';
	const char Constants::OTHER = 'O';
	const char Constants::WGS_CONTIG = 'W';
	
	// Gap types
	const std::string Constants::FRAGMENT = "fragment";
	const std::string Constants::CONTIG = "contig";
	const std::string Constants::SPLIT_FINISHED = "split_finished";
	const std::string Constants::CLONE = "clone";
	const std::string Constants::CENTROMERE = "centromere";
	const std::string Constants::SHORT_ARM = "short_arm";
	const std::string Constants::HETEROCHROMATIN = "heterochromatin";
	const std::string Constants::TELOMERE = "telomere";

	// Linkage strings
	const std::string Constants::LINKAGE_NO = "no";
	const std::string Constants::LINKAGE_YES = "yes";
	
} } }
