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

#ifndef __GFF_RECORD_HH__
#define __GFF_RECORD_HH__

#include <string>
#include <istream>
#include <ostream>
#include <list>

#include "bio/genome/BasicInterval.hh"
#include "util/io.hh"
#include "boost/unordered_map.hpp"

namespace bio { namespace gff {

	class GFFRecord {
	public:
		typedef float Score;

		struct Attribute {
			std::string tag;
			std::list<std::string> values;
		};
		
		GFFRecord();
		
		const std::string& getSeqname() const;
		const std::string& getSource() const;
		const std::string& getFeature() const;
		genome::Position getStart() const;
		genome::Position getEnd() const;
		Score getScore() const;
		genome::Strand getStrand() const;
		unsigned char getFrame() const;
		genome::BasicInterval getInterval() const;

		bool hasScore() const;
		bool hasStrand() const;
		bool hasFrame() const;

		void setInterval(const genome::Interval& i);
		void setSeqname(const std::string& seqname);
		void setSource(const std::string& source);
		void setFeature(const std::string& feature);
		void setStart(const genome::Position start);
		void setEnd(const genome::Position end);
		void setScore(const Score score);
		void setStrand(const genome::Strand strand);
		void setFrame(const unsigned char frame);

		void unsetScore();
		void unsetStrand();
		void unsetFrame();

		Attribute& getAttribute(const std::string& tag) const;
		bool hasAttribute(const std::string& tag) const;
		Attribute& addAttribute(const std::string& tag);
		void removeAttribute(const std::string& tag);

		friend std::ostream& operator<<(std::ostream& strm, const GFFRecord& r);
		friend std::istream& operator>>(std::istream& strm, GFFRecord& r);
		friend void operator>>(const std::string& tag, GFFRecord& r);

		bool operator<(const GFFRecord& r) const;
		
	private:
		typedef boost::unordered_map<std::string, const std::string*> StringMap;
		typedef std::list<Attribute> AttributeList;

		struct AttributeFinder {
			const std::string& tag;
			AttributeFinder(const std::string& tag);
			bool operator()(const Attribute& a) const;
		};

		const std::string* storeString(const std::string& s,
									   StringMap& strmap);
		void parseAttributes() const;
		
		static const std::string THE_EMPTY_STRING;
		static StringMap seqnames;
		static StringMap sources;
		static StringMap features;
		
		const std::string* seqname;
		const std::string* source;
		const std::string* feature;
		genome::Position start;
		genome::Position end;
		float score;
		char strand;
		char frame;
		mutable std::string attributeString;
		mutable AttributeList attributes;
		
		bool validScore;
	};

	inline GFFRecord::AttributeFinder::AttributeFinder(const std::string& tag)
		: tag(tag) {}
	inline bool GFFRecord::AttributeFinder::operator()(const Attribute& a) const {
		return a.tag == tag;
	}
	
	inline const std::string& GFFRecord::getSeqname() const { return *seqname; }
	inline const std::string& GFFRecord::getSource() const { return *source; }
	inline const std::string& GFFRecord::getFeature() const { return *feature; }
	inline genome::Position GFFRecord::getStart() const { return start; }
	inline genome::Position GFFRecord::getEnd() const { return end; }
	inline GFFRecord::Score GFFRecord::getScore() const { return score; }
	inline genome::Strand GFFRecord::getStrand() const { return strand; }
	inline unsigned char GFFRecord::getFrame() const { return frame - '0'; }

	inline bool GFFRecord::hasScore() const { return validScore; }
	inline bool GFFRecord::hasStrand() const { return strand != '.'; }
	inline bool GFFRecord::hasFrame() const { return frame != '.'; }
	
	inline void GFFRecord::setSeqname(const std::string& seqname) {
		this->seqname = storeString(seqname, seqnames);
	}
	inline void GFFRecord::setSource(const std::string& source) {
		this->source = storeString(source, sources);
	}
	inline void GFFRecord::setFeature(const std::string& feature) {
		this->feature = storeString(feature, features);
	}

	inline void GFFRecord::setStart(const genome::Position start) {
		this->start = start;
	}
	inline void GFFRecord::setEnd(const genome::Position end) {
		this->end = end;
	}

	inline void GFFRecord::setScore(const Score score) {
		this->score = score;
		this->validScore = true;
	}
	
	inline void GFFRecord::setStrand(const genome::Strand strand) {
		this->strand = strand;
	}
	inline void GFFRecord::setFrame(const unsigned char frame) {
		this->frame = '0' + frame;
	}

	inline void GFFRecord::unsetScore() { this->validScore = false; }
	inline void GFFRecord::unsetStrand() { this->strand = '.'; }
	inline void GFFRecord::unsetFrame() { this->frame = '.'; }

} }

#endif // __GFF_RECORD_HH__
