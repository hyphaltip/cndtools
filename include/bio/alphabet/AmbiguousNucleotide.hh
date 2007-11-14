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

#ifndef __AMBIGUOUS_NUCLEOTIDE_HH__
#define __AMBIGUOUS_NUCLEOTIDE_HH__

#include <string>
#include <vector>

#include "bio/alphabet/AmbiguousAlphabet.hh"
#include "bio/alphabet/Nucleotide.hh"

namespace bio { namespace alphabet {

	class AmbiguousNucleotide : public Nucleotide, public AmbiguousAlphabet {
	private:
		const std::vector<std::string> ambiguities;
		
	public:

		AmbiguousNucleotide(const std::string chars,
							const std::string complement,
							const std::vector<std::string> ambiguities);

		std::string getPossibleChars(const char c) const;
	};
	
	inline AmbiguousNucleotide::AmbiguousNucleotide(const std::string chars,
											 const std::string complement,
											 const std::vector<std::string> ambiguities)
		: Nucleotide(chars, complement), ambiguities(ambiguities) {
	}

	inline std::string AmbiguousNucleotide::getPossibleChars(const char c) const {
		return ambiguities[Nucleotide::encode(c)];
	}
	
	const std::string AMBIGUOUS_DNA_AMBIGUITIES[] =
		{"A", "C", "G", "T",
		 "ACGT",
		 "AC", "AG", "AT", "CG", "CT", "GT",
		 "ACG", "ACT", "AGT", "CGT"};

	const std::string AMBIGUOUS_RNA_AMBIGUITIES[] =
		{"A", "C", "G", "U",
		 "ACGU",
		 "AC", "AG", "AU", "CG", "CU", "GU",
		 "ACG", "ACU", "AGU", "CGU"};

	const AmbiguousNucleotide
	AmbiguousDNA("ACGTNMRWSYKVHDB",
				 "TGCANKYWSRMBDHV",
				 std::vector<std::string>(AMBIGUOUS_DNA_AMBIGUITIES,
										  AMBIGUOUS_DNA_AMBIGUITIES + sizeof(AMBIGUOUS_DNA_AMBIGUITIES)/sizeof(std::string)));
	
	const AmbiguousNucleotide
	GenomeDNA("ACGTN",
			  "TGCAN",
			  std::vector<std::string>(AMBIGUOUS_DNA_AMBIGUITIES,
									   AMBIGUOUS_DNA_AMBIGUITIES + 5));
	
	const AmbiguousNucleotide
	AmbiguousRNA("ACGUNMRWSYKVHDB",
				 "UGCANKYWSRMBDHV",
				 std::vector<std::string>(AMBIGUOUS_RNA_AMBIGUITIES,
										  AMBIGUOUS_RNA_AMBIGUITIES + sizeof(AMBIGUOUS_RNA_AMBIGUITIES)/sizeof(std::string)));

} }

#endif // __AMBIGUOUS_NUCLEOTIDE_HH__
