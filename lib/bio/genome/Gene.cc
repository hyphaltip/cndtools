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

#include "bio/genome/Gene.hh"

namespace bio { namespace genome {

	const Gene& Gene::getGene(const Exon& e) { return *e.gene; }
	const Gene& Gene::getGene(const Intron& i) { return *i.gene; }
	const Gene& Gene::getGene(const CDS& c) { return *c.gene; }
	const Gene& Gene::getGene(const UTR& u) { return *u.gene; }
	
	Intron Gene::get5PrimeIntron(const Exon& e) {
		if (e.gene->getStrand().isForward()) {
			return Intron(e.gene, e.i - 1);
		} else {
			return Intron(e.gene, e.i);
		}
	}
	
	Intron Gene::get3PrimeIntron(const Exon& e) {
		if (e.gene->getStrand().isForward()) {
			return Intron(e.gene, e.i);
		} else {
			return Intron(e.gene, e.i - 1);
		}
	}
	
	Exon Gene::get5PrimeExon(const Intron& i) {
		if (i.gene->getStrand().isForward()) {
			return Exon(i.gene, i.i);
		} else {
			return Exon(i.gene, i.i + 1);
		}
	}
	
	Exon Gene::get3PrimeExon(const Intron& i) {
		if (i.gene->getStrand().isForward()) {
			return Exon(i.gene, i.i + 1);
		} else {
			return Exon(i.gene, i.i);
		}
	}
	
	CDS Gene::getCDS(const Exon& e) {
		return CDS(e.gene, e.i);
	}
	
	Exon Gene::getExon(const CDS& c) {
		return Exon(c.gene, c.i);
	}
	
	Exon Gene::getExon(const UTR& u) {
		return Exon(u.gene, u.i);
	}

	UTR Gene::get5PrimeUTR(const Exon& e) {
		return UTR(e.gene, e.i, e.gene->getStrand().isForward());
	}

	UTR Gene::get3PrimeUTR(const Exon& e) {
		return UTR(e.gene, e.i, not e.gene->getStrand().isForward());
	}

	Gene::Gene() {}
	
	Gene::Gene(const std::string& name,
			   const std::string& chrom,
			   Strand strand,
			   Position codingStart,
			   Position codingEnd,
			   const std::vector<Position>& exonStarts,
			   const std::vector<Position>& exonEnds) :
		name(name),
		chrom(chrom),
		strand(strand),
		codingStart(codingStart),
		codingEnd(codingEnd),
		exonStarts(exonStarts),
		exonEnds(exonEnds),
		codingLength(0),
		codingStartIndex(exonStarts.size()),
		codingEndIndex(exonStarts.size())
	{
		initExons();
	}

	void Gene::set(const std::string& name,
				   const std::string& chrom,
				   Strand strand,
				   Position codingStart,
				   Position codingEnd,
				   const std::vector<Position>& exonStarts,
				   const std::vector<Position>& exonEnds) {
		this->name = name;
		this->chrom = chrom;
		this->strand = strand;
		this->codingStart = codingStart;
		this->codingEnd = codingEnd;
		this->exonStarts = exonStarts;
		this->exonEnds = exonEnds;
		this->codingLength = 0;
		this->codingStartIndex = exonStarts.size();
		this->codingEndIndex = exonEnds.size();

		initExons();
	}

	void Gene::initExons() {
		// Set indices of first and last coding exons
		for (size_t i = 0; i < exonStarts.size(); ++i) {
			if (exonStarts[i] < codingEnd and exonEnds[i] > codingStart) {
				if (codingStartIndex == exonStarts.size()) {
					codingStartIndex = i;
				}
				codingEndIndex = i;
			}
		}

		// Set phases and coding length
		phases.resize(getNumCDSs(), 0);
		for (size_t i = 0; i < phases.size(); ++i) {
			phases[i] = codingLength % 3;
			codingLength += getCDS(i).getLength();
		}
	}
	
	size_t Gene::getNumExons() const {
		return exonStarts.size();
	}

	size_t Gene::getNumIntrons() const {
		return getNumExons() - 1;
	}

	bool Gene::isCoding() const {
		return codingStartIndex != getNumExons();
	}
	
	size_t Gene::getNumCDSs() const {
		if (isCoding()) {
			return codingEndIndex - codingStartIndex + 1;
		} else {
			return 0;
		}
	}
	
	size_t Gene::getNumCDSIntrons() const {
		return getNumCDSs() - 1;
	}

	Exon Gene::getExon(size_t num) const {
		if (strand == '+') {
			return Exon(this, num);
		} else {
			return Exon(this, getNumExons() - 1 - num);
		}
	}
	
	CDS Gene::getCDS(size_t num) const {
		if (strand == '+') {
			return CDS(this, codingStartIndex + num);
		} else {
			return CDS(this, codingEndIndex - num);
		}
	}

	Intron Gene::getIntron(size_t num) const {
		if (strand == '+') {
			return Intron(this, num);
		} else {
			return Intron(this, getNumIntrons() - 1 - num);
		}
	}

	Position Gene::getTranscriptStart() const {
		return exonStarts.front();
	}
	
	Position Gene::getTranscriptEnd() const {
		return exonEnds.back();
	}
	
	Position Gene::getCodingStart() const {
		return codingStart;
	}
	
	Position Gene::getCodingEnd() const {
		return codingEnd;
	}
	
	Distance Gene::getTranscriptLength() const {
		return getTranscriptEnd() - getTranscriptStart();
	}
	
	Distance Gene::getCodingLength() const {
		return codingLength;
	}
	
	std::string Gene::getName() const {
		return name;
	}

	std::string Gene::getChrom() const {
		return chrom;
	}
	
	Position Gene::getStart() const {
		return getTranscriptStart();
	}
	
	Position Gene::getEnd() const {
		return getTranscriptEnd();
	}
	
	Strand Gene::getStrand() const {
		return strand;
	}

} }
