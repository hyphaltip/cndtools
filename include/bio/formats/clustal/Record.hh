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

#ifndef __BIO_FORMATS_CLUSTAL_RECORD_HH__
#define __BIO_FORMATS_CLUSTAL_RECORD_HH__

#include <string>
#include <vector>

#include "bio/alignment/NamedMultipleAlignment.hh"

namespace bio { namespace formats { namespace clustal {

    struct Record : public alignment::NamedMultipleAlignment {
		std::vector<std::string> names;
		std::vector<std::string> sequences;

		size_t getNumCols() const;
		size_t getNumSeqs() const;
		std::string getSeq(size_t seqNum) const;
		std::string getSubstring(size_t seqNum,
								 size_t colStart,
								 size_t len) const;
		char getChar(size_t seqNum, size_t colNum) const;

		std::string getName(size_t seqNum) const;
    };

} } }

#endif // __BIO_FORMATS_CLUSTAL_RECORD_HH__
