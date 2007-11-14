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

#ifndef __BIO_FORMATS_FASTA_OUTPUTSTREAM_HH__
#define __BIO_FORMATS_FASTA_OUTPUTSTREAM_HH__

#include <iosfwd>

#include "bio/formats/fasta/Record.hh"
#include "bio/formats/fasta/Constants.hh"
#include "bio/alignment/NamedMultipleAlignment.hh"

namespace bio { namespace formats { namespace fasta {

    class OutputStream : private Constants {
	public:
		OutputStream(std::ostream& strm);

		void setLineWidth(size_t lineWidth);

		OutputStream& operator<<(const Record& r);
		OutputStream& operator<<(const alignment::NamedMultipleAlignment& a);
		
	private:
		std::ostream& strm;
		size_t lineWidth;
    };

} } }

#endif // __BIO_FORMATS_FASTA_OUTPUTSTREAM_HH__
