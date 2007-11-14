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

#ifndef __BIO_FORMATS_NEWICK_HH__
#define __BIO_FORMATS_NEWICK_HH__

#include <iosfwd>

#include "bio/phylogenetic/Tree.hh"

namespace bio { namespace formats { namespace newick {

	// Read the first tree from STRM in newick format.  Returns NULL
	// if a tree can not be successfully read (either due to bad tree
	// syntax or IO error)
	phylogenetic::Tree* readTree(std::istream& strm);

	// Read the first tree in the string STR in newick format.  Returns NULL
	// if no tree can be read
	phylogenetic::Tree* readTree(const std::string& str);
	
	// Write the tree T onto STRM in newick format
	void writeTree(std::ostream& strm,
				   phylogenetic::Tree* t,
				   bool writeEdgeLengths = true);	

} } }

#endif // __BIO_FORMATS_NEWICK_HH__
