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

#ifndef __MATH_MAXPLUS_HH__
#define __MATH_MAXPLUS_HH__

#include <iosfwd>
#include <limits>
#include <algorithm>

namespace math {

	template<typename NumType>
	class MaxPlus {
	public:
		typedef MaxPlus Element;

		Element getZero() const;
		Element getMultiplicativeIdentity() const;
		
		MaxPlus(const NumType value = NumType());
	
		MaxPlus& operator*=(const MaxPlus& x);
		MaxPlus& operator+=(const MaxPlus& x);

		MaxPlus operator*(const MaxPlus& x) const;	
		MaxPlus operator+(const MaxPlus& x) const;

		bool operator<(const MaxPlus& x) const;
		bool operator<=(const MaxPlus& x) const;
		bool operator==(const MaxPlus& x) const;

		//		operator NumType() const;

		NumType value;
	};
		
	template<typename NumType>
	inline
	MaxPlus<NumType> MaxPlus<NumType>::getZero() const {
		return std::numeric_limits<NumType>::min();
	}

	template<>
	inline
	MaxPlus<double> MaxPlus<double>::getZero() const {
		return -std::numeric_limits<double>::max();
	}

	template<>
	inline
	MaxPlus<float> MaxPlus<float>::getZero() const {
		return -std::numeric_limits<float>::max();
	}
	
	template<typename NumType>
	inline
	MaxPlus<NumType> MaxPlus<NumType>::getMultiplicativeIdentity() const {
		return 0;
	}

	template<typename NumType>
	inline
	MaxPlus<NumType>::MaxPlus(const NumType value)
		: value(value) {
	}

	template<typename NumType>
	MaxPlus<NumType>&
	MaxPlus<NumType>::operator*=(const MaxPlus<NumType>& x) {
		if (*this == getZero() or x == getZero()) {
			*this = getZero();
		} else {
			value += x.value;
		}
		return *this;
	}

	template<typename NumType>
	inline MaxPlus<NumType>&
	MaxPlus<NumType>::operator+=(const MaxPlus<NumType>& x) {
		value = std::max(value, x.value);
		return *this;
	}

	template<typename NumType>
	inline MaxPlus<NumType>
	MaxPlus<NumType>::operator*(const MaxPlus<NumType>& x) const {
		return MaxPlus<NumType>(*this) *= x;
	}

	template<typename NumType>
	inline MaxPlus<NumType>
	MaxPlus<NumType>::operator+(const MaxPlus<NumType>& x) const {
		return MaxPlus<NumType>(*this) += x;
	}

	template<typename NumType>
	inline bool
	MaxPlus<NumType>::operator<(const MaxPlus<NumType>& x) const {
		return value < x.value;
	}

	template<typename NumType>
	inline bool
	MaxPlus<NumType>::operator==(const MaxPlus<NumType>& x) const {
		return value == x.value;
	}

	template<typename NumType>
	inline bool
	MaxPlus<NumType>::operator<=(const MaxPlus<NumType>& x) const {
		return value <= x.value;
	}
	
// 	template<typename NumType>
// 	inline
// 	MaxPlus<NumType>::operator NumType() const {
// 		return value;
// 	}

	template<typename NumType>
	inline std::ostream&
	operator<<(std::ostream& strm, const MaxPlus<NumType>& x) {
		return strm << x.value;
	}
			
}
		
#endif // __MATH_MAXPLUS_HH__
