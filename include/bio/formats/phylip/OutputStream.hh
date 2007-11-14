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

#ifndef __BIO_FORMATS_PHYLIP_OUTPUTSTREAM_HH__
#define __BIO_FORMATS_PHYLIP_OUTPUTSTREAM_HH__

#include <iosfwd>

#include "bio/formats/phylip/Constants.hh"
#include "bio/formats/phylip/Record.hh"

namespace bio { namespace formats { namespace phylip {

    class OutputStream : private Constants {
	public:
		OutputStream(std::ostream& strm);
		void setSequential(bool b);
		void setLineWidth(size_t lineWidth);
		void setSpaceFreq(size_t spaceFreq);
		void setBlockSpacing(size_t blockSpacing);
		OutputStream& operator<<(const Record& rec);
		operator bool() const;
		bool operator!() const;

	private:
		void writeBlankLines(size_t numLines);
		void writeFirstLine(const Record& rec);
		size_t writeSeq(const std::string& s,
						size_t colsWritten,
						size_t lineWidthRemaining);
		size_t writeSeqLine(const Record& rec,
							size_t i,
							size_t colsWritten);
		size_t writeFirstSeqLine(const Record& rec, size_t i);
		void writeInterleaved(const Record& rec);
		void writeSequential(const Record& rec);
		
		std::ostream& strm;
		bool hasWrittenRecord;
		bool sequential;
		size_t lineWidth;
		size_t spaceFreq;
		size_t blockSpacing;
	};

} } }

#endif // __BIO_FORMATS_PHYLIP_OUTPUTSTREAM_HH__
