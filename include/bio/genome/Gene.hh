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

#ifndef __BIO_GENOME_GENE_HH__
#define __BIO_GENOME_GENE_HH__

#include <vector>

#include "bio/genome/Interval.hh"
#include "bio/genome/Exon.hh"
#include "bio/genome/Intron.hh"
#include "bio/genome/CDS.hh"
#include "bio/genome/UTR.hh"

namespace bio { namespace genome {

    class Gene : public Interval {
	public:
		friend class Exon;
		friend class CDS;
		friend class Intron;
		friend class UTR;

		static const Gene& getGene(const Exon& e);
		static const Gene& getGene(const Intron& i);
		static const Gene& getGene(const CDS& c);
		static const Gene& getGene(const UTR& u);
		
		static Intron get5PrimeIntron(const Exon& e);
		static Intron get3PrimeIntron(const Exon& e);
		static Exon get5PrimeExon(const Intron& i);
		static Exon get3PrimeExon(const Intron& i);
		static CDS getCDS(const Exon& e);
		static Exon getExon(const CDS& c);
		static Exon getExon(const UTR& u);
		static UTR get5PrimeUTR(const Exon& e);
		static UTR get3PrimeUTR(const Exon& e);

		Gene();
		Gene(const std::string& name,
			 const std::string& chrom,
			 Strand strand,
			 Position codingStart,
			 Position codingEnd,
			 const std::vector<Position>& exonStarts,
			 const std::vector<Position>& exonEnds);

		void set(const std::string& name,
				 const std::string& chrom,
				 Strand strand,
				 Position codingStart,
				 Position codingEnd,
				 const std::vector<Position>& exonStarts,
				 const std::vector<Position>& exonEnds);
		
		size_t getNumExons() const;
		size_t getNumIntrons() const;
		size_t getNumCDSs() const;
		size_t getNumCDSIntrons() const;
		bool isCoding() const;
		Exon getExon(size_t i) const;
		CDS getCDS(size_t i) const;
		Intron getIntron(size_t i) const;
		Position getTranscriptStart() const;
		Position getTranscriptEnd() const;
		Position getCodingStart() const;
		Position getCodingEnd() const;
		Distance getTranscriptLength() const;
		Distance getCodingLength() const;

		std::string getName() const;
		
		std::string getChrom() const;
		Position getStart() const;
		Position getEnd() const;
		Strand getStrand() const;

	private:
		void initExons();

		std::string name;
		std::string chrom;
		Strand strand;
		Position codingStart;
		Position codingEnd;
		std::vector<Position> exonStarts;
		std::vector<Position> exonEnds;
		std::vector<int> phases;
		Distance codingLength;
		size_t codingStartIndex;
		size_t codingEndIndex;
    };

} }

#endif // __BIO_GENOME_GENE_HH__
