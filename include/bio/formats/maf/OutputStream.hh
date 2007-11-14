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

#ifndef __BIO_FORMATS_MAF_OUTPUTSTREAM_HH__
#define __BIO_FORMATS_MAF_OUTPUTSTREAM_HH__

#include <iosfwd>
#include <string>

#include "bio/formats/maf/Record.hh"
#include "bio/formats/maf/Header.hh"
#include "bio/formats/maf/Constants.hh"

namespace bio { namespace formats { namespace maf {

    class OutputStream : private Constants {
	public:
		OutputStream(std::ostream& strm, const Header& header = Header());
		OutputStream& operator<<(const Record& rec);
		OutputStream& operator<<(const std::string& s);
		operator bool() const;
		bool operator!() const;

	private:
		template<typename IntegerType>
		static size_t numDigits(IntegerType x);

		void writeHeader(const Header& header);
		void writeSequence(const Sequence& s,
						   size_t maxSrcLength,
						   size_t maxStartLength,
						   size_t maxSizeLength,
						   size_t maxSrcSizeLength);
		void writeVariables(const VariableMap& m);
		std::ostream& strm;
	};


	template<typename IntegerType>
	size_t
	OutputStream::numDigits(IntegerType x) {
		size_t digits = 0;
		while (x != 0) {
			++digits;
			x /= 10;
		}
		return digits;
	}
	
} } }


#endif // __BIO_FORMATS_MAF_OUTPUTSTREAM_HH__
