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

#ifndef __BIO_ALIGNMENT_BASICMULTIPLEALIGNMENT_HH__
#define __BIO_ALIGNMENT_BASICMULTIPLEALIGNMENT_HH__

#include <vector>
#include <string>

#include "bio/alignment/MultipleAlignment.hh"
#include "bio/alignment/Interval.hh"

namespace bio { namespace alignment {

    class BasicMultipleAlignment : public virtual MultipleAlignment {
	public:
		BasicMultipleAlignment();

		size_t getNumCols() const;
		size_t getNumSeqs() const;
		std::string getSeq(const size_t seqNum) const;

		std::string getSubstring(size_t seqNum,
								 size_t colStart,
								 size_t colEnd) const;
		char getChar(size_t seqNum, size_t colNum) const;

		void addSeq(const std::string& seq);
		void clear();
		
		size_t getColumnNum(const size_t seqNum,
							const size_t seqPos) const;
		Interval getColumnInterval(const size_t seqNum,
								   const size_t seqStart,
								   const size_t seqEnd) const;
		size_t getSeqPos(const size_t seqNum,
						 const size_t col) const;
		Interval getSeqInterval(const size_t seqNum,
								const size_t colStart,
								const size_t colEnd) const;
		std::string getCharsAfter(const size_t seqNum,
								  const size_t colStart,
								  const size_t numChars) const;
		std::string getCharsBefore(const size_t seqNum,
								   const size_t colEnd,
								   const size_t numChars) const;
		size_t getNumChars(const size_t seqNum,
						   const size_t colStart,
						   const size_t colEnd) const;

		void makeSeqsUppercase();
		
	private:
		std::vector<std::string> seqs;
	};

} }

#endif // __BIO_ALIGNMENT_BASICMULTIPLEALIGNMENT_HH__
