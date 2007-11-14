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

#ifndef __BIO_ALPHABET_NUCLEOTIDE_HH__
#define __BIO_ALPHABET_NUCLEOTIDE_HH__

#include <string>
#include <algorithm>

#include "util/charmanip.hh"
#include "util/stl.hh"
#include "bio/alphabet/Alphabet.hh"

namespace bio { namespace alphabet {

	class Nucleotide : public virtual Alphabet {
	private:
		static const util::charmanip::Mapper hardMasker;
		static const util::charmanip::Mapper unMasker;
		
		static const util::charmanip::Mapper transcriber;
		static const util::charmanip::Mapper reverseTranscriber;

		static const util::charmanip::Encoder nibEncoder;
		static const util::charmanip::Decoder nibDecoder;
		static const util::charmanip::Decoder nibHardMaskDecoder;
		static const util::charmanip::Decoder nibUnmaskDecoder;
		
		const util::charmanip::Decoder decoder;
		const util::charmanip::Encoder encoder;
		const util::charmanip::Mapper complementer;

	public:
		Nucleotide(const std::string chars, const std::string complement);
		
		char complement(const char base) const;
		
		unsigned char getSize() const;
		
		bool isMember(const char c) const;
		
		unsigned char encode(const char base) const;
		
		char decode(const unsigned char c) const;
		
		static std::string nibEncode(const std::string& seq);
		static std::string nibDecode(const std::string& seq,
									 bool hardMask=false,
									 bool unmask=false);

		void complementInPlace(std::string& seq) const;
		std::string complement(std::string seq) const;
		
		void reverseComplementInPlace(std::string& seq) const;
		std::string reverseComplement(std::string seq) const;

		static void transcribeInPlace(std::string& seq);
	    static std::string transcribe(std::string seq);

		static void reverseTranscribeInPlace(std::string& seq);
	    static std::string reverseTranscribe(std::string seq);

		static void hardMaskInPlace(std::string& seq);
		static std::string hardMask(std::string seq);

		static void unMaskInPlace(std::string& seq);
		static std::string unMask(std::string seq);
	};

	inline Nucleotide::Nucleotide(const std::string chars, const std::string complement) :
		decoder(chars),
		encoder(chars, false),
		complementer(chars, complement, false)
	{}

	inline char Nucleotide::complement(const char base) const {
		return complementer(base);
	}
		
	inline unsigned char Nucleotide::getSize() const {
		return decoder.decoding.size();
	}
		
	inline bool Nucleotide::isMember(const char c) const {
		return encoder.isMember(c);
	}
		
	inline unsigned char Nucleotide::encode(const char base) const {
		return encoder(base);
	}
		
	inline char Nucleotide::decode(const unsigned char c) const {
		return decoder(c);
	}
	
	//const Nucleotide AmbiguousDNA("ACGTNMRWSYKVHDB", "TGCANKYWSRMBDHV");
	//const Nucleotide GenomeDNA("ACGTN", "TGCAN");
	const Nucleotide DNA("ACGT", "TGCA");
	const Nucleotide DNA_PURINES("AG", "TC");
	const Nucleotide DNA_PYRIMIDINES("CT", "GA");
	//const Nucleotide AmbiguousRNA("ACGUNMRWSYKVHDB", "UGCANKYWSRMBDHV");
	const Nucleotide RNA("ACGU", "UGCA");
	const Nucleotide RNA_PURINES("AG", "UC");
	const Nucleotide RNA_PYRIMIDINES("CU", "GA");

} }

#endif // __BIO_ALPHABET_NUCLEOTIDE_HH__
