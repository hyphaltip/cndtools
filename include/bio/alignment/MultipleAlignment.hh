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

#ifndef __BIO_ALIGNMENT_MULTIPLEALIGNMENT_HH__
#define __BIO_ALIGNMENT_MULTIPLEALIGNMENT_HH__

#include <string>
#include <vector>

namespace bio { namespace alignment {

	class MultipleAlignment {
	public:
		virtual size_t getNumCols() const = 0;
		virtual size_t getNumSeqs() const = 0;
		virtual std::string getSeq(size_t seqNum) const = 0;

		virtual std::string getSubstring(size_t seqNum,
										 size_t colStart,
										 size_t colEnd) const;
		
		virtual char getChar(size_t seqNum, size_t colNum) const;

		virtual size_t getSeqLen(const size_t seqNum) const;
		
		virtual size_t getSeqLen(const size_t seqNum,
								 const size_t colStart,
								 const size_t colEnd) const;

		virtual bool isColIdentical(const size_t colNum) const;

		virtual void getGapMask(const size_t colNum,
								std::vector<bool>& mask) const;
		
		virtual size_t getNumBlocks() const;
		
		virtual ~MultipleAlignment();
	};

} }

#endif // __BIO_ALIGNMENT_MULTIPLEALIGNMENT_HH__
