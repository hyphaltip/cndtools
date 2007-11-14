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

#ifndef __BIO_ALIGNMENT_NAMEDMULTIPLEALIGNMENT_HH__
#define __BIO_ALIGNMENT_NAMEDMULTIPLEALIGNMENT_HH__

#include <string>

#include "bio/alignment/MultipleAlignment.hh"

namespace bio { namespace alignment {

    class NamedMultipleAlignment : public virtual MultipleAlignment {
	public:
		// Returns the name of the sequence at index I
		virtual std::string getName(size_t i) const = 0;

		// Return the index of the sequence with name NAME or -1 if
		// there is no such sequence
		virtual int getSeqNum(const std::string& name) const;
    };

} }

#endif // __BIO_ALIGNMENT_NAMEDMULTIPLEALIGNMENT_HH__
