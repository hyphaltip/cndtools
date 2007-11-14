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

#include "VectorMultiSet.hh"

#include <ostream_wrapper.h>
using namespace polymake;

VectorMultiSet& VectorMultiSet::operator+=(const VectorMultiSet& s) {
	EltCountList unionEltCounts;
	EltCountList::const_iterator pos1 = eltCounts.begin();
	EltCountList::const_iterator pos2 = s.eltCounts.begin();

	// Process the elements of both sets in increasing order
	while (pos1 != eltCounts.end() and pos2 != s.eltCounts.end()) {
		if (pos1->first < pos2->first) {
			unionEltCounts.push_back(*pos1);
			++pos1;
		} else if (pos2->first < pos1->first) {
			unionEltCounts.push_back(*pos2);
			++pos2;
		} else {
			unionEltCounts.push_back(EltCount(pos1->first,
											  pos1->second + pos2->second));
			++pos1;
			++pos2;
		}
	}

	// Add remaining vectors from set that we did not get to the end of
	if (pos1 != eltCounts.end()) {
		EltCountList::const_iterator end = eltCounts.end();
		unionEltCounts.insert(unionEltCounts.end(), pos1, end);
	} else if (pos2 != s.eltCounts.end()) {
		unionEltCounts.insert(unionEltCounts.end(), pos2, s.eltCounts.end());
	}

	eltCounts.swap(unionEltCounts);

	return *this;
}

VectorMultiSet& VectorMultiSet::operator*=(const VectorMultiSet& s) {
	// Modify elements in place if one of the sets has only one element
	if (eltCounts.size() == 1 or s.eltCounts.size() == 1) {
		EltCount shift;
		if (s.eltCounts.size() > 1) {
			shift = eltCounts.front();
			eltCounts = s.eltCounts;
		} else {
			shift = s.eltCounts.front();
		}
		for (EltCountList::iterator pos = eltCounts.begin();
			 pos != eltCounts.end(); ++pos) {
			pos->first += shift.first;
			pos->second *= shift.second;
		}
	}
	// Both sets have more than one element, therefore we need to
	// add points pairwise, sort, and then combine identical points
	else {
		// Add points pairwise
		EltCountList productEltCounts;
		for (EltCountList::const_iterator pos1 = eltCounts.begin();
			 pos1 != eltCounts.end(); ++pos1) {
			for (EltCountList::const_iterator pos2 = s.eltCounts.begin();
				 pos2 != s.eltCounts.end(); ++pos2) {
				EltCount product(pos1->first + pos2->first,
								 pos1->second * pos2->second);
				productEltCounts.push_back(product);
			}
		}

		// Sort all points
		std::sort(productEltCounts.begin(), productEltCounts.end());

		// Combine identical points and put into this multiset
		eltCounts.clear();
		for (EltCountList::const_iterator pos = productEltCounts.begin();
			 pos != productEltCounts.end(); ++pos) {
			if (eltCounts.empty() or pos->first != eltCounts.back().first) {
				eltCounts.push_back(*pos);
			} else {
				eltCounts.back().second += pos->second;
			}
		}
			
	}
	return *this;
}

std::ostream& operator<<(std::ostream& strm,
						 const VectorMultiSet& s) {
	for (VectorMultiSet::const_iterator it = s.begin(); it != s.end(); ++it) {
		wrap(strm) << it->first;
		strm << ' ' << it->second << '\n';
	}
	return strm;
}
