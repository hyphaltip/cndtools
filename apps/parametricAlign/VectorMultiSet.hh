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

#ifndef __VECTOR_MULTI_SET_HH__
#define __VECTOR_MULTI_SET_HH__

#include <vector>
#include <ostream>

#include <Vector.h>
#include <Rational.h>
using namespace polymake;

class VectorMultiSet {
public:
	typedef Vector<Rational> Element;
	typedef size_t Count;
	typedef std::pair<Element, Count> EltCount;
	typedef std::vector<EltCount> EltCountList;
	typedef EltCountList::const_iterator const_iterator;

	// Initialize an empty multiset
	VectorMultiSet();

	// Initialize a multiset with one element ELT
	VectorMultiSet(const Element& elt);

	static VectorMultiSet makeVariable(size_t var_index, size_t dims);
	static VectorMultiSet makeConstant(Rational val, size_t dims);
	
	const_iterator begin() const;
	const_iterator end() const;

	// Add VectorMultiSet S to this multiset.  This multiset will now
	// be the union of the two multisets.
	VectorMultiSet& operator+=(const VectorMultiSet& s);

	// Multiply this multiset by the VectorMultiSet S.  This multiset will now
	// be the minkowski sum of the two multisets.
	VectorMultiSet& operator*=(const VectorMultiSet& s);
	
	VectorMultiSet operator*(const VectorMultiSet& s) const;
	VectorMultiSet operator+(const VectorMultiSet& s) const;

	friend std::ostream& operator<<(std::ostream& strm,
									const VectorMultiSet& s);

private:
	EltCountList eltCounts;
};


// Initialize an empty multiset
inline VectorMultiSet::VectorMultiSet() {}

// Initialize a multiset with one element ELT
inline VectorMultiSet::VectorMultiSet(const Element& elt) {
	eltCounts.push_back(EltCount(elt, 1));
}

inline VectorMultiSet::const_iterator VectorMultiSet::begin() const {
	return eltCounts.begin();
}

inline VectorMultiSet::const_iterator VectorMultiSet::end() const {
	return eltCounts.end();
}
	
inline VectorMultiSet VectorMultiSet::operator*(const VectorMultiSet& s) const {
	VectorMultiSet result = *this;
	return (result *= s);
}

inline VectorMultiSet VectorMultiSet::operator+(const VectorMultiSet& s) const {
	VectorMultiSet result = *this;
	return (result += s);
}

inline VectorMultiSet VectorMultiSet::makeVariable(size_t var_index,
												   size_t dims) {
	VectorMultiSet::Element v(dims);
	v[var_index] = 1;
	return VectorMultiSet(v);
}

inline VectorMultiSet VectorMultiSet::makeConstant(Rational val, size_t dims) {
	if (val == 0) {
		return VectorMultiSet();
	} else {
		VectorMultiSet::Element c(dims);
		return VectorMultiSet(c);
	}
}


#endif // __VECTOR_MULTI_SET_HH__
