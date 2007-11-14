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

#ifndef __BIO_ALPHABET_AMINOACID_HH__
#define __BIO_ALPHABET_AMINOACID_HH__

#include <string>

#include "bio/alphabet/Alphabet.hh"
#include "util/charmanip.hh"

namespace bio { namespace alphabet {

	class AminoAcid : public Alphabet {
	private:
		static const util::charmanip::Mapper hardMasker;
		static const util::charmanip::Mapper unMasker;
		
		const util::charmanip::Decoder decoder;
		const util::charmanip::Encoder encoder;

	public:
		AminoAcid(const std::string chars);
		
		unsigned char getSize() const;
		
		bool isMember(const char c) const;
		
		unsigned char encode(const char base) const;
		
		char decode(const unsigned char c) const;
		
		static void hardMaskInPlace(std::string& seq);
		static std::string hardMask(std::string seq);

		static void unMaskInPlace(std::string& seq);
		static std::string unMask(std::string seq);
	};

	const AminoAcid AmbiguousAminoAcid("ACDEFGHIKLMNPQRSTVWYBXZ");
	const AminoAcid UnambiguousAminoAcid("ACDEFGHIKLMNPQRSTVWY");
	
} }

#endif // __BIO_ALPHABET_AMINOACID_HH__
