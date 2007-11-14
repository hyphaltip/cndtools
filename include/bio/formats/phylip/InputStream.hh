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

#ifndef __BIO_FORMATS_PHYLIP_INPUTSTREAM_HH__
#define __BIO_FORMATS_PHYLIP_INPUTSTREAM_HH__

#include <iosfwd>
#include <string>

#include "util/io/line/InputStream.hh"
#include "bio/formats/phylip/Constants.hh"
#include "bio/formats/phylip/Record.hh"

namespace bio { namespace formats { namespace phylip {

	class InputStream : private Constants {
	public:
		InputStream(std::istream& strm);
		void setSequential(bool sequential);

		InputStream& operator>>(Record& rec);
		operator bool() const;
		bool operator!() const;

	private:
		void initRec(Record& rec, size_t numSeqs, size_t numCols) const;
		void appendSequence(std::string::const_iterator start,
							std::string::const_iterator end,
							std::string& s) const;
		void readUntilNextNonBlankLine();
		void readFirstLine(size_t& numRecs, size_t& numCols);
		void readFirstSeqLine(Record& rec, size_t i);
		void readSeqLine(Record& rec, size_t i);
		void readInterleaved(Record& rec, size_t numCols);
		void readSequential(Record& rec, size_t numCols);

		util::io::line::InputStream strm;
		bool sequential;
		std::string line;
	};

} } }

#endif // __BIO_FORMATS_PHYLIP_INPUTSTREAM_HH__
