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
#include <cstdlib> // For atoi, atol, atof
#include "boost/tokenizer.hpp"

#include "bio/gff/GFFRecord.hh"
#include "util/parser.hh"

namespace bio { namespace gff {

	GFFRecord::GFFRecord() :
		seqname(&THE_EMPTY_STRING),
		source(&THE_EMPTY_STRING),
		feature(&THE_EMPTY_STRING),
		start(0),
		end(0),
		score(0),
		strand('.'),
		frame('.'),
		attributeString(""),
		attributes(),
		validScore(false)
	{}

	const std::string GFFRecord::THE_EMPTY_STRING = "";
	GFFRecord::StringMap GFFRecord::seqnames = GFFRecord::StringMap();
	GFFRecord::StringMap GFFRecord::sources = GFFRecord::StringMap();
	GFFRecord::StringMap GFFRecord::features = GFFRecord::StringMap();

	bool GFFRecord::operator<(const GFFRecord& r) const {
		int chromComp = getSeqname().compare(r.getSeqname());
		return chromComp < 0
			or (chromComp == 0
				and (getStart() < r.getStart()
					 or (getStart() == r.getStart()
						 and getEnd() < r.getEnd())));
	}
	
	void GFFRecord::parseAttributes() const {
		// Check if attribute parsing has already been done (or if
		// there are no attributes to parse)
		if (attributeString.empty()) {
			return;
		}
		
		typedef boost::tokenizer<boost::char_separator<char> > CharTok;
		typedef boost::tokenizer<boost::escaped_list_separator<char> > EscListTok;

		CharTok ctok(attributeString, boost::char_separator<char>(";"));
		for (CharTok::iterator pos = ctok.begin(); pos != ctok.end(); ++pos) {
			EscListTok ltok(*pos,
							boost::escaped_list_separator<char>('\\', ' ', '\"'));
			Attribute a;
			for (EscListTok::iterator sub = ltok.begin();
				 sub != ltok.end(); ++sub) {
				if (sub->empty()) {
					continue;
				} else if (a.tag.empty()) {
					a.tag = *sub;
				} else {
					a.values.push_back(*sub);
				}
			}
			if (!a.tag.empty()) {
				attributes.push_back(a);
			}
		}

		// Make the attribute string empty to free memory and indicate
		// that attribute parsing has been done
		attributeString.clear();
	}

	GFFRecord::Attribute& GFFRecord::getAttribute(const std::string& tag) const {
		parseAttributes();
		AttributeList::iterator pos = std::find_if(attributes.begin(),
												   attributes.end(),
												   AttributeFinder(tag));
		if (pos == attributes.end()) {
			throw std::runtime_error(std::string("Attribute not found: ")
									 + tag);
		} else {
			return *pos;
		}
	}

	bool GFFRecord::hasAttribute(const std::string& tag) const {
		parseAttributes();
		return std::find_if(attributes.begin(),
							attributes.end(),
							AttributeFinder(tag)) != attributes.end();
	}

	GFFRecord::Attribute& GFFRecord::addAttribute(const std::string& tag) {
		parseAttributes();
		attributes.push_back(Attribute());
		attributes.back().tag = tag;
		return attributes.back();
	}

	void GFFRecord::removeAttribute(const std::string& tag) {
		parseAttributes();
		attributes.erase(std::remove_if(attributes.begin(), attributes.end(),
										AttributeFinder(tag)),
						 attributes.end());
	}

	const std::string*
	GFFRecord::storeString(const std::string& s, StringMap& strmap) {
		StringMap::const_iterator pos = strmap.find(s);
		if (pos == strmap.end()) {
			pos = strmap.insert(std::make_pair(s, new std::string(s))).first;
		}
		return pos->second;
	}

	std::ostream& operator<<(std::ostream& strm, const GFFRecord& r) {

		strm << r.getSeqname() << '\t'
			 << r.getSource() << '\t'
			 << r.getFeature() << '\t'
			 << r.getStart() << '\t'
			 << r.getEnd() << '\t';

		if (r.hasScore()) {
			strm << r.getScore() << '\t';
		} else {
			strm << '.' << '\t';
		}

		strm << r.strand << '\t'
			 << r.frame << '\t';

		if (r.attributeString.empty()) {
			GFFRecord::AttributeList::const_iterator pos;
			for (pos = r.attributes.begin(); pos != r.attributes.end(); ++pos) {
				if (pos != r.attributes.begin()) {
					strm << " ; ";
				}
				strm << pos->tag;
				std::list<std::string>::const_iterator vpos;
				for (vpos = pos->values.begin();
					 vpos != pos->values.end(); ++vpos) {
					if (util::parser::isFreeText(*vpos)) {
						strm << ' ' <<  '"' << *vpos << '"';
					} else {
						strm << ' ' << *vpos;
					}
				}
			}
		} else {
			strm << r.attributeString;
		}

		strm << '\n';

		return strm;
	}

	std::istream& operator>>(std::istream& strm, GFFRecord& rec) {
		// Ignore whitespace and comments
		while (true) {
			util::parser::ignoreWhitespace(strm);
			
			// If we have hit EOF, return
			if (strm.eof()) {
				strm.ignore();
				return strm;
			}
			
			// If we have hit a comment line, skip rest of line and continue
			// otherwise, break because we are at a real record line
			if (strm.peek() == '#') {
				util::parser::ignoreLine(strm);
			} else {
				break;
			}
		}

		// Read in the entire record line
		std::string line;
		if (std::getline(strm, line)) {
			line >> rec;
		}

		return strm;
	}

	void operator>>(const std::string& line, GFFRecord& rec) {
		// Split the string at tabs
		typedef boost::char_separator<char> CharSep;
		typedef boost::tokenizer<CharSep> CharTok;
		CharTok tabtok(line, CharSep("\t", "", boost::keep_empty_tokens));

		// Step through tokens
		CharTok::iterator t = tabtok.begin();

		// Parse seqname
 		if (t == tabtok.end()) {
 			throw util::parser::FormatError(line, "invalid number of fields");
 		}
		rec.setSeqname(*t);

		// Parse source
		++t;
		if (t == tabtok.end()) {
			throw util::parser::FormatError(line, "invalid number of fields");
		}
		rec.setSource(*t);

		// Parse feature
		++t;
		if (t == tabtok.end()) {
			throw util::parser::FormatError(line, "invalid number of fields");
		}		
		rec.setFeature(*t);

		// Parse start
		++t;
		if (t == tabtok.end()) {
			throw util::parser::FormatError(line, "invalid number of fields");
		}
		rec.start = std::atoi(t->c_str());

		// Parse end
		++t;
		if (t == tabtok.end()) {
			throw util::parser::FormatError(line, "invalid number of fields");
		}
		rec.end = std::atoi(t->c_str());

		// Parse score
		++t;
		if (t == tabtok.end()) {
			throw util::parser::FormatError(line, "invalid number of fields");
		}		
		if (*t == ".") {
			rec.unsetScore();
		} else {
			rec.setScore(std::atof(t->c_str()));
		}

		// Parse strand
		++t;
		if (t == tabtok.end()) {
			throw util::parser::FormatError(line, "invalid number of fields");
		}
		rec.strand = (*t).at(0);

		// Parse frame
		++t;
		if (t == tabtok.end()) {
			throw util::parser::FormatError(line, "invalid number of fields");
		}
		rec.frame = (*t).at(0);
		
		// Parse attributes
		++t;
		if (t == tabtok.end()) {
			throw util::parser::FormatError(line, "invalid number of fields");
		}
		std::string attributeString = *t;
		++t;
		while (t != tabtok.end()) {
			attributeString += " " + *t;
			++t;
		}
		// Remove previous attributes
		rec.attributes.clear();
		// Strip off any comments
		std::string::size_type commentBegin = attributeString.find('#');
		if (commentBegin != std::string::npos) {
			rec.attributeString = attributeString.substr(0, commentBegin);
		} else {
			rec.attributeString = attributeString;
		}

	}

	genome::BasicInterval GFFRecord::getInterval() const {
		return genome::BasicInterval(*seqname,
									 start - 1,
									 end,
									 (hasStrand() ? strand : '+'));
	}

	void GFFRecord::setInterval(const genome::Interval& i) {
		setSeqname(i.getChrom());
		setStart(i.getStart() + 1);
		setEnd(i.getEnd());
		setStrand(i.getStrand());
	}

} }
