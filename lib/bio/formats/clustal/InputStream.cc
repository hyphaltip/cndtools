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
#include <istream>

#include "bio/formats/clustal/InputStream.hh"
#include "util/string.hh"
using util::string::startsWith;

namespace bio { namespace formats { namespace clustal {

	InputStream::InputStream(std::istream& strm)
		: strm(strm),
		  hasNewRecord(false) {
	}
	
	InputStream::operator bool() const { return hasNewRecord or strm; }

	bool InputStream::operator!() const { return not *this; }

	void InputStream::readFirstLine() {
		strm >> line;
		if (not strm or not startsWith(line, HEADER_LINE_PREFIX)) {
			throw std::runtime_error("Invalid first line:\n" + line);
		}
	}

	void InputStream::initRec(Record& rec) const {
		rec.names.clear();
		rec.sequences.clear();
	}

	void InputStream::readBlocks(Record& rec) {
		bool firstBlock = true;
		std::vector<std::string> lines;
		while (true) {
			lines.clear();
			readBlockLines(lines);
			if (lines.empty()) { break; }
			parseBlockLines(lines, rec, firstBlock);
			firstBlock = false;
		}
	}

	void InputStream::readBlockLines(std::vector<std::string>& lines) {
		while (strm >> line and line.empty()) {}
		if (not strm) { return; }
		do {
			lines.push_back(line);
		} while (strm >> line and not line.empty());
	}

	void InputStream::parseBlockLine(const std::string& line,
									 std::string& name,
									 std::string& seq) const {
		std::istringstream lineStream(line);
		lineStream >> name >> seq;
		if (not lineStream) {
			throw std::runtime_error("Invalid alignment line:\n" + line);
		}
	}
	
	void InputStream::parseBlockLines(std::vector<std::string>& lines,
									  Record& rec,
									  bool first) const {
		if (not first and lines.size() - 1 != rec.getNumSeqs()) {
			throw std::runtime_error("Block does not have same number of lines "
									 "as first block");
		}
		std::string name, seq;
		for (size_t i = 0; i < lines.size() - 1; ++i) {
			parseBlockLine(lines[i], name, seq);
			if (first) {
				rec.names.push_back(name);
				rec.sequences.push_back(seq);
			} else {
				rec.sequences[i].append(seq);
			}
		}
	}
	
	InputStream& InputStream::operator>>(Record& rec) {
		if (hasNewRecord) {
			hasNewRecord = false;
			return *this;
		}
		readFirstLine();
		initRec(rec);
		readBlocks(rec);
		hasNewRecord = not rec.names.empty();
		return *this;
	}

} } }
