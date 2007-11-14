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

#ifndef __BIO_FORMATS_CLUSTAL_OUTPUTSTREAM_HH__
#define __BIO_FORMATS_CLUSTAL_OUTPUTSTREAM_HH__

#include <vector>
#include <iosfwd>

#include "bio/alignment/NamedMultipleAlignment.hh"
#include "bio/formats/clustal/Constants.hh"
#include "bio/formats/clustal/Record.hh"

namespace bio { namespace formats { namespace clustal {

    class OutputStream : private Constants {
	public:
		OutputStream(std::ostream& strm);
		void setLineWidth(size_t lineWidth);
		void setMinNameSeqSpacing(size_t minNameSeqSpacing);
		void setMinGutterLen(size_t minGutterLen);
		void setBlockSpacing(size_t blockSpacing);
		void setSeqNos(bool b);
		OutputStream& operator<<(const Record& rec);
		OutputStream& operator<<(const alignment::NamedMultipleAlignment& align);
		operator bool() const;
		bool operator!() const;

	private:
		void writeBlankLines(size_t numLines);
		void writeFirstLine();
		void writeBlockLine(const std::string& name,
							const std::string& seq,
							size_t gutterLen,
							size_t& seqPos);
		void writeBlock(const Record& rec,
						size_t start,
						size_t gutterLen,
						size_t seqLen,
						std::vector<size_t>& seqPos);
		void writeIdentityLine(const std::vector<std::string>& v,
							   size_t start,
							   size_t seqLen,
							   size_t gutterLen);
		void writeBlocks(const Record& rec);
		
		static bool isIdenticalColumn(const std::vector<std::string>& v,
									  size_t pos);
		static size_t maxLen(const std::vector<std::string>& v);
		static size_t maxSeqNosLen(const std::vector<std::string>& v);
		static size_t numDigits(size_t n);
		
		std::ostream& strm;
		size_t lineWidth;
		size_t minNameSeqSpacing;
		size_t minGutterLen;
		size_t blockSpacing;
		bool seqnos;
	};

} } }

#endif // __BIO_FORMATS_CLUSTAL_OUTPUTSTREAM_HH__
