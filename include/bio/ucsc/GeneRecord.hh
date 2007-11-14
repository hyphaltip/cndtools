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

#ifndef __GENE_RECORD_HH__
#define __GENE_RECORD_HH__

#include <string>
#include <istream>
#include <ostream>
#include <vector>
#include <functional>

#include "bio/genome/Interval.hh"
#include "util/string.hh"

namespace bio { namespace ucsc {

	class GeneRecord {
	public:
		GeneRecord();

		template<class Iterator>
		GeneRecord(const std::string& name,
				   const std::string& chrom,
				   const genome::Strand strand,
				   const genome::Position txStart,
				   const genome::Position txEnd,
				   const genome::Position cdsStart,
				   const genome::Position cdsEnd,
				   const size_t exonCount,
				   Iterator exonStartsBegin, Iterator exonStartsEnd,
				   Iterator exonEndsBegin, Iterator exonEndsEnd);

		const std::string& getName() const;
		const std::string& getChrom() const;
		genome::Strand getStrand() const;
		genome::Position getTxStart() const;
		genome::Position getTxEnd() const;
		genome::Position getCDSStart() const;
		genome::Position getCDSEnd() const;
		size_t getExonCount() const;
		genome::Position getExonStart(const size_t exonNum) const;
		genome::Position getExonEnd(const size_t exonNum) const;
		
		friend std::ostream& operator<<(std::ostream& strm, const GeneRecord& r);
		friend std::istream& operator>>(std::istream& strm, GeneRecord& r);

	private:
		std::string name;
		std::string chrom;
		genome::Strand strand;
		genome::Position txStart;
		genome::Position txEnd;
		genome::Position cdsStart;
		genome::Position cdsEnd;
		size_t exonCount;
		std::vector<genome::Position> exonStarts;
		std::vector<genome::Position> exonEnds;
	};

	inline GeneRecord::GeneRecord()
		: name(),
		  chrom(),
		  strand('+'),
		  txStart(0),
		  txEnd(0),
		  cdsStart(0),
		  cdsEnd(0),
		  exonCount(0),
		  exonStarts(),
		  exonEnds() {
	}

	template<class Iterator>
	GeneRecord::GeneRecord(const std::string& name,
						   const std::string& chrom,
						   const genome::Strand strand,
						   const genome::Position txStart,
						   const genome::Position txEnd,
						   const genome::Position cdsStart,
						   const genome::Position cdsEnd,
						   const size_t exonCount,
						   Iterator exonStartsBegin, Iterator exonStartsEnd,
						   Iterator exonEndsBegin, Iterator exonEndsEnd) 
		: name(name),
		  chrom(chrom),
		  strand(strand),
		  txStart(txStart),
		  txEnd(txEnd),
		  cdsStart(cdsStart),
		  cdsEnd(cdsEnd),
		  exonCount(exonCount),
		  exonStarts(exonStartsBegin, exonStartsEnd),
		  exonEnds(exonEndsBegin, exonEndsEnd) {
	}
		
	inline const std::string& GeneRecord::getName() const { return name; }

	inline const std::string& GeneRecord::getChrom() const { return chrom; }

	inline genome::Strand GeneRecord::getStrand() const { return strand; }

	inline genome::Position GeneRecord::getTxStart() const {
		return txStart;
	}

	inline genome::Position GeneRecord::getTxEnd() const {
		return txEnd;
	}

	inline genome::Position GeneRecord::getCDSStart() const {
		return cdsStart;
	}

	inline genome::Position GeneRecord::getCDSEnd() const {
		return cdsEnd;
	}

	inline size_t GeneRecord::getExonCount() const { return exonCount; }

	inline genome::Position
	GeneRecord::getExonStart(const size_t exonNum) const {
		return exonStarts[exonNum];
	}

	inline genome::Position
	GeneRecord::getExonEnd(const size_t exonNum) const {
		return exonEnds[exonNum];
	}

	inline std::ostream& operator<<(std::ostream& strm, const GeneRecord& r) {
		strm << r.name << '\t'
			 << r.chrom << '\t'
			 << r.strand << '\t'
			 << r.txStart << '\t'
			 << r.txEnd << '\t'
			 << r.cdsStart << '\t'
			 << r.cdsEnd << '\t'
			 << r.exonCount << '\t';

		std::vector<genome::Position>::const_iterator pos;
		for (pos = r.exonStarts.begin(); pos != r.exonStarts.end(); ++pos) {
			strm << *pos << ',';
		}
		strm << '\t';
		for (pos = r.exonEnds.begin(); pos != r.exonEnds.end(); ++pos) {
			strm << *pos << ',';
		}
			
		return strm << '\n';
	}
	
	inline std::istream& operator>>(std::istream& strm, GeneRecord& r) {
		const size_t MIN_NUM_FIELDS = 10;
		std::string line;
		if (!std::getline(strm, line)) {
			return strm;
		}
		std::vector<std::string> fields;
		fields.reserve(MIN_NUM_FIELDS);
		util::string::split(line, std::back_inserter(fields), "\t");
		if (fields.size() < MIN_NUM_FIELDS) {
			throw std::runtime_error("Incorrect # of fields: " + line);
		}

		util::string::Converter<genome::Position> toDistance;
		util::string::Converter<size_t> toSizeT;
			
		r.name = fields[0];
		r.chrom = fields[1];
		r.strand = fields[2].at(0);
		r.txStart = toDistance(fields[3]);
		r.txEnd = toDistance(fields[4]);
		r.cdsStart = toDistance(fields[5]);
		r.cdsEnd = toDistance(fields[6]);
		r.exonCount = toSizeT(fields[7]);

		r.exonStarts.clear();
		r.exonEnds.clear();
		r.exonStarts.reserve(r.exonCount);
		r.exonEnds.reserve(r.exonCount);
			
		typedef boost::tokenizer<boost::char_separator<char> > CharTok;
			
		CharTok starts(fields[8], boost::char_separator<char>(","));
		std::transform(starts.begin(), starts.end(),
					   std::back_inserter(r.exonStarts),
					   toDistance);
		CharTok ends(fields[9], boost::char_separator<char>(","));
		std::transform(ends.begin(), ends.end(),
					   std::back_inserter(r.exonEnds),
					   toDistance);
			
		return strm;
	}

} }

#endif // __GENE_RECORD_HH__
