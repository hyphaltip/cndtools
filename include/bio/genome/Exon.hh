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

#ifndef __BIO_GENOME_EXON_HH__
#define __BIO_GENOME_EXON_HH__

#include "bio/genome/Interval.hh"

namespace bio { namespace genome {

	class Gene;
	
    class Exon : public Interval {
	public:
		friend class Gene;
		
		enum Type {
			SINGLE,
			INITIAL,
			INTERNAL,
			TERMINAL,
			NUM_TYPES
		};
		
		std::string getChrom() const;
		Position getStart() const;
		Position getEnd() const;
		Strand getStrand() const;

		size_t getNum() const;
		
		bool hasCDS() const;
		bool has5PrimeUTR() const;
		bool has3PrimeUTR() const;
		bool has5PrimeIntron() const;
		bool has3PrimeIntron() const;

		Type getType() const;
		
	private:
		Exon(const Gene* gene, size_t i);

		const Gene* gene;
		size_t i;
    };

} }

#endif // __BIO_GENOME_EXON_HH__
