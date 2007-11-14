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

#include "bio/genome/Exon.hh"
#include "bio/genome/Gene.hh"

namespace bio { namespace genome {

	std::string Exon::getChrom() const {
		return gene->getChrom();
	}
	
	Position Exon::getStart() const {
		return gene->exonStarts[i];
	}
	
	Position Exon::getEnd() const {
		return gene->exonEnds[i];
	}
	
	Strand Exon::getStrand() const {
		return gene->getStrand();
	}
	
	size_t Exon::getNum() const {
		return (getStrand() == '+' ? i : gene->getNumExons() - 1 - i);
	}

	bool Exon::hasCDS() const {
		return (i >= gene->codingStartIndex and i <= gene->codingEndIndex);
	}
		
	bool Exon::has5PrimeUTR() const {
		return (getStrand() == '+'
				? getStart() < gene->codingStart
				: getEnd() > gene->codingEnd);
	}
	
	bool Exon::has3PrimeUTR() const {
		return (getStrand() == '+'
				? getEnd() > gene->codingEnd
				: getStart() < gene->codingStart);
	}
	
	bool Exon::has5PrimeIntron() const {
		return i != 0;
	}
	
	bool Exon::has3PrimeIntron() const {
		return i != (gene->getNumExons() - 1);
	}

	Exon::Type Exon::getType() const {
		if (gene->getNumExons() == 1) {
			return SINGLE;
		} else {
			const size_t num = getNum();
			if (num == 0) {
				return INITIAL;
			} else if (num == (gene->getNumExons() - 1)) {
				return TERMINAL;
			} else {
				return INTERNAL;
			}
		}
	}
		
	Exon::Exon(const Gene* gene, size_t i) : gene(gene), i(i) {
	}
	
} }
