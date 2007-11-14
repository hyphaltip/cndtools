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

#include <sstream>
#include <algorithm>

#include "bio/ucsc/GeneInputStream.hh"
#include "util/string.hh"

namespace bio { namespace ucsc {

	const size_t MIN_NUM_FIELDS = 10;
	
	GeneInputStream::GeneInputStream(std::istream& strm,
									 const size_t bufferSize)
		: strm(strm, bufferSize) {
	}
	
	GeneInputStream::operator bool() { return strm; }
	bool GeneInputStream::operator!() { return !strm; }

	GeneInputStream& GeneInputStream::operator>>(genome::Gene& g) {
		strm >> line;
		if (not strm) { return *this; }

		std::string name;
		std::string chrom;
		genome::Strand strand;
		genome::Position txStart, txEnd, cdsStart, cdsEnd;
		size_t exonCount;
		std::string exonStartsString, exonEndsString;
		
		std::istringstream stringStream(line);
		stringStream >> name >> chrom >> strand
					 >> txStart >> txEnd >> cdsStart >> cdsEnd
					 >> exonCount >> exonStartsString >> exonEndsString;
		if (not stringStream) { 
			throw std::runtime_error("Invalid line:\n" + line);
		}

		std::vector<genome::Position> exonStarts, exonEnds;
		exonStarts.reserve(exonCount);
		exonEnds.reserve(exonCount);

		typedef boost::tokenizer<boost::char_separator<char> > CharTok;
		util::string::Converter<genome::Position> toPosition;
			
		CharTok starts(exonStartsString, boost::char_separator<char>(","));
		std::transform(starts.begin(), starts.end(),
					   std::back_inserter(exonStarts),
					   toPosition);
		CharTok ends(exonEndsString, boost::char_separator<char>(","));
		std::transform(ends.begin(), ends.end(),
					   std::back_inserter(exonEnds),
					   toPosition);

		g.set(name, chrom, strand, cdsStart, cdsEnd, exonStarts, exonEnds);
		
		return *this;
	}
	
} }
