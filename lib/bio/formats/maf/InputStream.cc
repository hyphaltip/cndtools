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

#include <sstream>

#include "bio/formats/maf/InputStream.hh"
#include "bio/formats/maf/Record.hh"
#include "util/string.hh"
using util::string::startsWith;
using util::string::split;

namespace bio { namespace formats { namespace maf {

	InputStream::InputStream(std::istream& strm)
		: strm(strm), hasNewRecord(false) {
		readHeader();
	}
	
	InputStream::operator bool() const { return strm or hasNewRecord; }
	bool InputStream::operator!() const { return not (strm or hasNewRecord); }

	const Header& InputStream::getHeader() const {
		return header;
	}

	void InputStream::parseVariables(const std::string& s, VariableMap& m) {
		std::vector<std::string> tokens;
		split(s, std::back_inserter(tokens));
		typedef std::vector<std::string>::const_iterator Iterator;
		for (Iterator i = tokens.begin(); i != tokens.end(); ++i) {
			const std::string& token = *i;
			size_t splitPoint = token.find('=');
			if (splitPoint == std::string::npos) {
				throw std::runtime_error("Invalid variable/value pair: "
										 + token);
			}
			m.setVariable(token.substr(0, splitPoint),
						  token.substr(splitPoint + 1));
		}
	}
	
	void InputStream::readHeader() {
		strm >> line;
		if (not startsWith(line, HEADER_LINE_PREFIX)) {
			throw std::runtime_error("Invalid MAF header line:\n" + line);
		}
		parseVariables(line.substr(HEADER_LINE_PREFIX.size()), header);
		if (not header.hasVariable("version")) {
			throw std::runtime_error("Invalid MAF header line, "
									 "version not defined:\n" + line);
		}
		strm >> line;
		if (not startsWith(line, COMMENT_LINE_PREFIX)) {
			throw std::runtime_error("Invalid MAF parameter line:\n" + line);
		}
		header.parameters = line.substr(COMMENT_LINE_PREFIX.size());
	}

	void InputStream::skipToNextParagraph() {
		while (strm >> line) {
			if (not (line.empty() or startsWith(line, COMMENT_LINE_PREFIX))) {
				break;
			}
		}
	}		
	
	InputStream& InputStream::operator>>(Record& rec) {
		hasNewRecord = false;
		// Read paragraphs until an alignment block is found
		do {
			skipToNextParagraph();
		} while (strm and not startsWith(line, ALIGNMENT_LINE_PREFIX));
		if (not strm) { return *this; }
		
		// Parse alignment line
		rec.variables.clear();
		parseVariables(line.substr(ALIGNMENT_LINE_PREFIX.size()), rec);

		// Parse sequence lines
		rec.sequences.clear();
		while (strm >> line) {
			// Stop when end of paragraph has been reached
			if (line.empty()) { break; }

			// Skip over non-sequence lines
			if (not startsWith(line, SEQ_LINE_PREFIX)) {
				continue;
			}

			// Create a new sequence from the line
			rec.sequences.push_back(Sequence());
			Sequence& seq = rec.sequences.back();
			std::istringstream seqStream(line.substr(SEQ_LINE_PREFIX.size()));
			seqStream >> seq.src >> seq.start >> seq.size >> seq.strand
					  >> seq.srcSize >> seq.text;
			if (not seqStream) {
				throw std::runtime_error("Invalid MAF sequence line:\n" + line);
			}
		}
		hasNewRecord = true;
		return *this;
	}
	
} } }
