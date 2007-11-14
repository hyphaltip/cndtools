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

#include "bio/genome/CDS.hh"
#include "bio/genome/Gene.hh"

namespace bio { namespace genome {

	std::string CDS::getChrom() const {
		return gene->getChrom();
	}
	
	Position CDS::getStart() const {
		if (i > gene->codingStartIndex) {
			return gene->exonStarts[i];
		} else {
			return gene->codingStart;
		}
	}
	
	Position CDS::getEnd() const {
		if (i < gene->codingEndIndex) {
			return gene->exonEnds[i];
		} else {
			return gene->codingEnd;
		}
	}
	
	Strand CDS::getStrand() const {
		return gene->getStrand();
	}

	size_t CDS::getNum() const {
		if (getStrand() == '+') {
			return i - gene->codingStartIndex;
		} else {
			return gene->codingEndIndex - i;
		}
	}
	
	int CDS::getPhase() const {
		return gene->phases[getNum()];
	}

	CDS::CDS(const Gene* gene, size_t i) : gene(gene), i(i) {
	}

} }
