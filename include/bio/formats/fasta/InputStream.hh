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

#ifndef __BIO_FORMATS_FASTA_INPUTSTREAM_HH__
#define __BIO_FORMATS_FASTA_INPUTSTREAM_HH__

#define FASTA_BUFFER_SIZE (4 * 1024)

#include <iosfwd>
#include <string>
#include <vector>

#include "bio/formats/fasta/Constants.hh"
#include "bio/formats/fasta/Record.hh"
#include "bio/alignment/BasicNamedMultipleAlignment.hh"

namespace bio { namespace formats { namespace fasta {

	class InputStream : private Constants {
	public:
		InputStream(std::istream& strm,
					const size_t bufferSize=FASTA_BUFFER_SIZE);
		
		// If given a value of TRUE, whitespace (including newlines)
		// are kept inside of the sequence part of each FASTA record.
		// If given a value of FALSE, whitespace is removed from the
		// sequence portion of each record.
		void setKeepWhitespace(bool b);

		InputStream& operator>>(Record& rec);
		InputStream& operator>>(alignment::BasicNamedMultipleAlignment& a);
		operator bool();
		bool operator!();
 	
		typedef Record ValueType;

	private:
		void fillBuffer();		

		std::istream& strm;
		std::vector<char> buffer;
		std::vector<char>::iterator pos;
		bool atEnd;
		bool keepWhitespace;
	};


} } }

#endif // __BIO_FORMATS_FASTA_INPUTSTREAM_HH__
