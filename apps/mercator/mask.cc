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

#include "mask.hh"

unsigned int firstInMask(const Mask& m) {
	unsigned int first;
	for (first = 0; first < m.size() && !m.test(first); ++first);
	return first;
}

void makeSubsetMasks(Mask m, vector<Mask>& subsets,
					 unsigned int minSubsetSize) {
	unsigned int count = m.count();
	if (count == minSubsetSize) {
		subsets.push_back(m);
	} else if (count > minSubsetSize) {
		subsets.reserve(1 << count);
		unsigned int first = firstInMask(m);
		m.reset(first);
		vector<Mask> withoutFirst;
		makeSubsetMasks(m, withoutFirst,
						minSubsetSize ? minSubsetSize - 1 : 0);
		vector<Mask>::const_iterator s;
		for (s = withoutFirst.begin(); s != withoutFirst.end(); ++s) {
			if (s->count() >= minSubsetSize) {
				subsets.push_back(*s);
			}
			subsets.push_back(*s);
			subsets.back().set(first);
		}
	}
}

bool MaskSorter::operator()(const Mask& m1, const Mask& m2) const {
	return (m1.count() < m2.count() ||
			(m1.count() == m2.count() &&
			 m1.to_ulong() < m2.to_ulong()));
}

void makeMasks(const int numGenomes, vector<Mask>& masks) {
	makeSubsetMasks(Mask((1 << numGenomes) - 1),
					masks, 2);
	sort(masks.rbegin(), masks.rend(), MaskSorter());
}
