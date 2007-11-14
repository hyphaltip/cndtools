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

#include "bio/alphabet/Nucleotide.hh"
#include "bio/alphabet/Alphabet.hh"
#include "util/charmanip.hh"
#include "util/stl.hh"

namespace bio {

namespace alphabet {
		
	std::string Nucleotide::nibEncode(const std::string& seq) {
		unsigned int length = seq.size() % 2 ?
			seq.size() / 2 + 1 : seq.size() / 2;
		std::string encoded(length, 0);
		for (unsigned int i = 0; i < seq.size(); i += 2) {
			unsigned char upper = nibEncoder(seq[i]);
			unsigned char lower = (i + 1 < seq.size()) ?
				nibEncoder(seq[i + 1]) : 0x0F;
			encoded[i / 2] = (upper << 4) | lower;
		}
		return encoded;
	}

	std::string Nucleotide::nibDecode(const std::string& seq,
									  bool hardMask,
									  bool unmask) {
		const util::charmanip::Decoder* nd = &nibDecoder;
		if (hardMask) {
			nd = &nibHardMaskDecoder;
		} else if (unmask) {
			nd = &nibUnmaskDecoder;
		}
									 
		unsigned int length;
		if (seq.size() > 0 && ((seq[seq.size() - 1] & 0x0F) == 0x0F)) {
			length = seq.size() * 2 - 1;
		} else {
			length = seq.size() * 2;
		}
		std::string decoded(length, '?');
		for (unsigned int i = 0; i < length; i += 2) {
			decoded[i] = (*nd)((seq[i / 2] >> 4) & 0x0F);
		}
		for (unsigned int i = 1; i < length; i += 2) {
			decoded[i] = (*nd)(seq[i / 2] & 0x0F);
		}
		return decoded;
	}

	void Nucleotide::complementInPlace(std::string& seq) const {
		std::transform(seq.begin(), seq.end(), seq.begin(),
					   util::stl::make_functor_ref(complementer));
	}
		
	std::string Nucleotide::complement(std::string seq) const {
		complementInPlace(seq);
		return seq;
	}
		
	void Nucleotide::reverseComplementInPlace(std::string& seq) const {
		std::reverse(seq.begin(), seq.end());
		complementInPlace(seq);
	}
		
	std::string Nucleotide::reverseComplement(std::string seq) const {
		reverseComplementInPlace(seq);
		return seq;
	}

	void Nucleotide::transcribeInPlace(std::string& seq) {
		std::transform(seq.begin(), seq.end(), seq.begin(),
					   util::stl::make_functor_ref(transcriber));
	}

	std::string Nucleotide::transcribe(std::string seq) {
		transcribeInPlace(seq);
		return seq;
	}

	void Nucleotide::reverseTranscribeInPlace(std::string& seq) {
		std::transform(seq.begin(), seq.end(), seq.begin(),
					   util::stl::make_functor_ref(reverseTranscriber));
	}

	std::string Nucleotide::reverseTranscribe(std::string seq) {
		reverseTranscribeInPlace(seq);
		return seq;
	}

	void Nucleotide::hardMaskInPlace(std::string& seq) {
		std::transform(seq.begin(), seq.end(), seq.begin(),
					   util::stl::make_functor_ref(hardMasker));
	}

	std::string Nucleotide::hardMask(std::string seq) {
		hardMaskInPlace(seq);
		return seq;
	}

	void Nucleotide::unMaskInPlace(std::string& seq) {
		std::transform(seq.begin(), seq.end(), seq.begin(),
					   util::stl::make_functor_ref(unMasker));
	}

	std::string Nucleotide::unMask(std::string seq) {
		unMaskInPlace(seq);
		return seq;
	}

	const util::charmanip::Mapper
	Nucleotide::hardMasker("abcdefghijklmnopqrstuvwxyz",
						   "NNNNNNNNNNNNNNNNNNNNNNNNNN");
	const util::charmanip::Mapper
	Nucleotide::unMasker("abcdefghijklmnopqrstuvwxyz",
						 "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	
	const util::charmanip::Mapper
	Nucleotide::transcriber("T", "U", false);
	const util::charmanip::Mapper
	Nucleotide::reverseTranscriber("U", "T", false);
	
	const util::charmanip::Encoder Nucleotide::nibEncoder("ACGTUNacgtun-", true, 'N');
	const util::charmanip::Decoder Nucleotide::nibDecoder("ACGTUNacgtun-");
	const util::charmanip::Decoder Nucleotide::nibHardMaskDecoder("ACGTUNNNNNNN-");
	const util::charmanip::Decoder Nucleotide::nibUnmaskDecoder("ACGTUNACGTUN-");
};
};
