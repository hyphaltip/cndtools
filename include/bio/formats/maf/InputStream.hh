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

#ifndef __BIO_FORMATS_MAF_INPUTSTREAM_HH__
#define __BIO_FORMATS_MAF_INPUTSTREAM_HH__

#include <iosfwd>

#include "bio/formats/maf/types.hh"
#include "bio/formats/maf/Header.hh"
#include "bio/formats/maf/Constants.hh"
#include "util/io/line/InputStream.hh"

namespace bio { namespace formats { namespace maf {

	class InputStream : private Constants {
	public:
		InputStream(std::istream& strm);
		
		const Header& getHeader() const;

		InputStream& operator>>(Record& rec);
		operator bool() const;
		bool operator!() const;

	private:
		void readHeader();
		void skipToNextParagraph();
		static void parseVariables(const std::string& s, VariableMap& m);
		
		util::io::line::InputStream strm;
		std::string line;
		Header header;
		bool hasNewRecord;
	};

} } }

#endif // __BIO_FORMATS_MAF_INPUTSTREAM_HH__
