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

#ifndef __MASK_HH__
#define __MASK_HH__

#include "types.hh"

// Returns the index of the first bit set in mask M.
// If no bit is set, returns the size of M
unsigned int firstInMask(const Mask& m);

void makeSubsetMasks(Mask m, vector<Mask>& subsets,
					 unsigned int minSubsetSize = 0);

void makeMasks(const int numGenomes, vector<Mask>& masks);

struct MaskSorter {
	bool operator()(const Mask& m1, const Mask& m2) const;
};

#endif // __MASK_HH__
