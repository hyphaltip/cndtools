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

#ifndef __BIO_FORMATS_AGP_RECORD_HH__
#define __BIO_FORMATS_AGP_RECORD_HH__

#include <string>
#include <iosfwd>
#include <vector>
#include "boost/lexical_cast.hpp"

#include "bio/genome/BasicInterval.hh"
#include "util/io.hh"
#include "util/string.hh"

// Facilities for processing AGP formatted records
// See http://genome.ucsc.edu/goldenPath/datorg.html for information

namespace bio { namespace formats { namespace agp {

	class Record {
	public:

		Record();

		// Source sequence constructor
		Record(const std::string& chrom,
			   const std::string& contig,
			   const genome::Position contigStart,
			   const genome::Position contigEnd,
			   const unsigned int agpRecNum,
			   const char sourceType,
			   const std::string& sourceAccession,
			   const genome::Position sourceStart,
			   const genome::Position sourceEnd,
			   const char orientation);
		
		// Gap constructor
		Record(const std::string& chrom,
			   const std::string& contig,
			   const genome::Position contigStart,
			   const genome::Position contigEnd,
			   const unsigned int agpRecNum,			   
			   const char sourceType,
			   const genome::Distance gapLength,
			   const std::string& gapKind,
			   const bool gapBridged);

		// Return the chromosome that is being assembled from source
		// sequences
		std::string getChrom() const;

		// Return the contig that is being assembled from source
		// sequences.  Returns the empty string if only a chromosome
		// is specified
		std::string getContig() const;

		// Return the start position of where the source sequence is
		// placed in the chromosome/contig
		genome::Position getContigStart() const;

		// Return the end position of where the source sequence is
		// placed in the chromosome/contig
		genome::Position getContigEnd() const;

		// Return the interval of the source sequence
		genome::BasicInterval getTargetInterval() const;
		
		// Return the length of the interval in the contig
		genome::Distance getContigLength() const;
		
		// Returns the AGP record number
		unsigned int getRecNum() const;

		// Return the type of sequence that this record represents One
		// of: F: Finished, A: in Active finishing, D: Draft, P:
		// PreDraft, O: Other, N: gap
		char getSourceType() const;

		// Returns true if this is a gap record (type == 'N')
		bool isGap() const;

		// Set the interval of the target
		void setTargetInterval(const genome::Interval& i);
		
		// Set the chromosome that is being assembled from source
		// sequences
		void setChrom(const std::string& chrom);

		// Set the contig that is being assembled from source
		// sequences.  Set as the empty string if only a chromosome
		// is specified
		void setContig(const std::string& contig);

		// Set the start position of where the source sequence is
		// placed in the chromosome/contig
		void setContigStart(genome::Position start);

		// Set the end position of where the source sequence is
		// placed in the chromosome/contig
		void setContigEnd(genome::Position end);

		// Set the AGP record number
		void setRecNum(unsigned int num);

		// Set the type of sequence that this record represents One
		// of: F: Finished, A: in Active finishing, D: Draft, P:
		// PreDraft, O: Other, N: gap
		void setSourceType(char type);

		// Set this record as a gap record (type == 'N')
		void setAsGap();

		
		// Methods only available for source (non-gap) records

		// Return the interval of the source sequence
		genome::BasicInterval getSourceInterval() const;
		
		// Returns the accession.version string for the source sequence
		std::string getSourceAccession() const;

		// Returns the start in the source sequence
		genome::Position getSourceStart() const;

		// Returns the end in the source sequence
		genome::Position getSourceEnd() const;

		// Returns the length of the interval in the source sequence
		genome::Distance getSourceLength() const;
		
		// Returns the orientation ('+' or '-') of the source sequence
		// in the contig
		char getSourceOrientation() const;

		
		// Methods only available for gap records

		// Returns the length of the gap (number of Ns to be placed in
		// the chromosome/contig)
		genome::Distance getGapLength() const;

		// Returns the kind of gap:
		//
		// fragment - a gap between two sequence contigs (also called
		// a "sequence gap")
		//
		// split_finished - a special sized gap between two finished
		// sequence contigs
		//
		// clone - a gap between two clones that do not overlap
		//
		// contig - a gap between clone contigs in the genome layout
		// (also called a "layout gap")
		//
		// centromere - a gap inserted for the centromere
		//
		// short_arm - a gap inserted at the start of an acrocentric
		// chromosome
		//
		// heterochromatin - a gap inserted for an especially large
		// region of heterochromatin (may include the centromere as
		// well.)
		//
		// telomere - a gap inserted for a telomere
		std::string getGapKind() const;

		// Returns true if if there is a cDNA or BACend pair or
		// plasmid end pair that spans the gap
		bool isGapBridged() const;

		// Set the interval of the source
		void setSourceInterval(const genome::Interval& i);
		void setSourceAccession(const std::string& accession);
		void setSourceStart(genome::Position start);
		void setSourceEnd(genome::Position end);
		void setSourceOrientation(char orientation);
		
		void setGapLength(genome::Distance length);
		void setGapKind(const std::string& kind);
		void setGapBridged(bool bridged);

		friend std::ostream& operator<<(std::ostream& strm, const Record& r);
		friend std::istream& operator>>(std::istream& strm, Record& r);
		friend void operator>>(const std::string& line, Record& r);

		static void stripComment(std::string& line);

	private:
		std::vector<std::string> fields;
	};

	inline Record::Record() {
		fields.reserve(9); // There can be at most 9 fields
	}
	
	inline std::string Record::getChrom() const {
		std::string::size_type pos = fields[0].find('/');
		return (pos == std::string::npos ? fields[0] : fields[0].substr(0, pos));
	}
		
	inline std::string Record::getContig() const {
		std::string::size_type pos = fields[0].find('/');
		return (pos == std::string::npos ? "" : fields[0].substr(pos));
	}
	
	inline genome::Position Record::getContigStart() const {
		return boost::lexical_cast<genome::Position>(fields[1]);
	}
	
	inline genome::Position Record::getContigEnd() const {
		return boost::lexical_cast<genome::Position>(fields[2]);
	}

	inline genome::Distance Record::getContigLength() const {
		return getContigEnd() - getContigStart() + 1;
	}

	inline genome::BasicInterval Record::getTargetInterval() const {
		return genome::BasicInterval(getChrom(),
									 getContigStart() - 1,
									 getContigEnd());
	}
	
	inline unsigned int Record::getRecNum() const {
		return boost::lexical_cast<unsigned int>(fields[3]);
	}

	inline char Record::getSourceType() const {
		return fields[4][0];
	}
	
	inline bool Record::isGap() const {
		return getSourceType() == 'N';
	}

	inline genome::BasicInterval Record::getSourceInterval() const {
		assert(!isGap());
		return genome::BasicInterval(getSourceAccession(),
									 getSourceStart() - 1,
									 getSourceEnd(),
									 getSourceOrientation());
	}
	
	inline std::string Record::getSourceAccession() const {
		assert(!isGap());
		return fields[5];
	}
		
	inline genome::Position Record::getSourceStart() const {
		assert(!isGap());
		return boost::lexical_cast<genome::Position>(fields[6]);
	}
	
	inline genome::Position Record::getSourceEnd() const {
		assert(!isGap());
		return boost::lexical_cast<genome::Position>(fields[7]);
	}

	inline genome::Distance Record::getSourceLength() const {
		return getSourceEnd() - getSourceStart() + 1;
	}
	
	inline char Record::getSourceOrientation() const {
		assert(!isGap());
		return fields[8][0];
	}
		
	inline genome::Distance Record::getGapLength() const {
		assert(isGap());
		return boost::lexical_cast<genome::Distance>(fields[5]);
	}
		
	inline std::string Record::getGapKind() const {
		assert(isGap());
		return fields[6];
	}
	
	inline bool Record::isGapBridged() const {
		assert(isGap());
		return fields[7] == "yes";
	}

	inline void Record::setTargetInterval(const genome::Interval& i) {
		setChrom(i.getChrom());
		setContigStart(i.getStart());
		setContigEnd(i.getEnd());
	}
	
	inline void Record::setChrom(const std::string& chrom) {
		std::string::size_type pos = fields[0].find('/');
		if (pos == std::string::npos) {
			fields[0] = chrom;
		} else {
			fields[0] = chrom + fields[0].substr(pos);
		}
	}
	
	inline void Record::setContig(const std::string& contig) {
		std::string::size_type pos = fields[0].find('/');
		if (pos == std::string::npos) {
			fields[0] += "/" + contig;
		} else {
			fields[0] = fields[0].substr(0, pos + 1) + contig;
		}
	}		
	
	inline void Record::setContigStart(genome::Position start) {
		fields[1] = util::string::toString(start);
	}

	inline void Record::setContigEnd(genome::Position end) {
		fields[2] = util::string::toString(end);
	}		
	
	inline void Record::setRecNum(unsigned int num) {
		fields[3] = util::string::toString(num);
	}		

	inline void Record::setSourceType(char type) {
		fields[4] = util::string::toString(type);
		fields.resize(isGap() ? 8 : 9);
	}
	
	inline void Record::setAsGap() {
		setSourceType('N');
	}

	inline void Record::setSourceInterval(const genome::Interval& i) {
		setSourceAccession(i.getChrom());
		setSourceStart(i.getStart());
		setSourceEnd(i.getEnd());
	}
	
	inline void Record::setSourceAccession(const std::string& accession) {
		assert(!isGap());
		fields[5] = accession;
	}
		
	inline void Record::setSourceStart(genome::Position start) {
		assert(!isGap());
		fields[6] = util::string::toString(start);
	}
	
	inline void Record::setSourceEnd(genome::Position end) {
		assert(!isGap());
		fields[7] = util::string::toString(end);
	}
	
	inline void Record::setSourceOrientation(char orientation) {
		assert(!isGap());
		fields[8] = util::string::toString(orientation);
	}
		

	inline void Record::setGapLength(genome::Distance length) {
		assert(isGap());
		fields[5] = util::string::toString(length);
	}
		
	inline void Record::setGapKind(const std::string& kind) {
		assert(isGap());
		fields[6] = kind;
	}
	
	inline void Record::setGapBridged(bool bridged) {
		assert(isGap());
		fields[7] == (bridged ? "yes" : "no");
	}

} } }

#endif // __BIO_FORMATS_AGP_RECORD_HH__
