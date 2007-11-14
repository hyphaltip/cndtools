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

#ifndef __BIO_FORMATS_CLUSTAL_INPUTSTREAM_HH__
#define __BIO_FORMATS_CLUSTAL_INPUTSTREAM_HH__

#include <iosfwd>
#include <string>

#include "util/io/line/InputStream.hh"
#include "bio/formats/clustal/Constants.hh"
#include "bio/formats/clustal/Record.hh"

namespace bio { namespace formats { namespace clustal {

	class InputStream : private Constants {
	public:
		InputStream(std::istream& strm);

		InputStream& operator>>(Record& rec);
		operator bool() const;
		bool operator!() const;

	private:
		void readFirstLine();
		void initRec(Record& rec) const;
		void readBlocks(Record& rec);
		void readBlockLines(std::vector<std::string>& lines);
		void parseBlockLine(const std::string& line,
										 std::string& name,
										 std::string& seq) const;
		void parseBlockLines(std::vector<std::string>& lines,
							 Record& rec,
							 bool first) const;

		util::io::line::InputStream strm;
		std::string line;
		bool hasNewRecord;
	};

} } }

#endif // __BIO_FORMATS_CLUSTAL_INPUTSTREAM_HH__
