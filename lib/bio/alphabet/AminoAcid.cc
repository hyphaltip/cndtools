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

#include <string>
#include <algorithm>

#include "bio/alphabet/AminoAcid.hh"
#include "bio/alphabet/Alphabet.hh"
#include "util/charmanip.hh"
#include "util/stl.hh"

namespace bio {

namespace alphabet {

	AminoAcid::AminoAcid(const std::string chars) 
		: decoder(chars),
		  encoder(chars, false)
	{}

	unsigned char AminoAcid::getSize() const {
		return decoder.decoding.size();
	}
		
	bool AminoAcid::isMember(const char c) const {
		return encoder.isMember(c);
	}
		
	unsigned char AminoAcid::encode(const char base) const {
		return encoder(base);
	}
		
	char AminoAcid::decode(const unsigned char c) const {
		return decoder(c);
	}
		
	void AminoAcid::hardMaskInPlace(std::string& seq) {
		std::transform(seq.begin(), seq.end(), seq.begin(),
					   util::stl::make_functor_ref(hardMasker));
	}

	std::string AminoAcid::hardMask(std::string seq) {
		hardMaskInPlace(seq);
		return seq;
	}

	void AminoAcid::unMaskInPlace(std::string& seq) {
		std::transform(seq.begin(), seq.end(), seq.begin(),
					   util::stl::make_functor_ref(unMasker));
	}

	std::string AminoAcid::unMask(std::string seq) {
		unMaskInPlace(seq);
		return seq;
	}

	const util::charmanip::Mapper
	AminoAcid::hardMasker("abcdefghijklmnopqrstuvwxyz",
						  "XXXXXXXXXXXXXXXXXXXXXXXXXX");
	const util::charmanip::Mapper
	AminoAcid::unMasker("abcdefghijklmnopqrstuvwxyz",
						"ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	
};

};
