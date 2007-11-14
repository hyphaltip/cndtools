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

#ifndef __BIO_GENOME_STRAND_HH__
#define __BIO_GENOME_STRAND_HH__

#include <string>
#include <ostream>
#include <istream>

namespace bio { namespace genome {

	class Strand {
	public:
		Strand(const char strand='+');
		void flip();
		Strand opposite() const;
		operator char() const;
		bool isForward() const;
		bool isReverse() const;

		friend std::ostream& operator<<(std::ostream& strm, const Strand& s);
		friend std::istream& operator>>(std::istream& strm, Strand& s);
	private:
		char strand;
	};

	inline Strand::Strand(const char strand) : strand(strand) {}
	inline void Strand::flip() { *this = opposite(); }
	inline Strand Strand::opposite() const {
		return (strand == '+' ? '-' : '+');
	}
	inline Strand::operator char() const { return strand; }
	inline bool Strand::isForward() const { return strand == '+'; }
	inline bool Strand::isReverse() const { return strand == '-'; }
	inline std::ostream& operator<<(std::ostream& strm, const Strand& s) {
		return strm << s.strand;
	}
	inline std::istream& operator>>(std::istream& strm, Strand& s) {
		std::string strandString;
		strm >> strandString;
		if (strm) { s.strand = strandString.at(0); }
		return strm;
	}
	
} }

#endif // __BIO_GENOME_STRAND_HH__
