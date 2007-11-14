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

#ifndef __MAX_TIMES_HH__
#define __MAX_TIMES_HH__

#include <algorithm>

namespace math {

	template<typename NumType>
	class MaxTimes {
	public:
		typedef MaxTimes<NumType> Element;
		static const NumType zero;
		static const NumType multiplicativeIdentity;
	
		NumType value;

		MaxTimes<NumType>(const NumType value = NumType());
	
		MaxTimes<NumType>& operator*=(const MaxTimes<NumType>& x);
		MaxTimes<NumType>& operator+=(const MaxTimes<NumType>& x);

		MaxTimes<NumType> operator*(const MaxTimes<NumType>& x) const;	
		MaxTimes<NumType> operator+(const MaxTimes<NumType>& x) const;

		bool operator<(const MaxTimes<NumType>& x) const;
		bool operator==(const MaxTimes<NumType>& x) const;
	};

	template<typename NumType>
	const NumType MaxTimes<NumType>::zero = 0;

	template<typename NumType>
	const NumType MaxTimes<NumType>::multiplicativeIdentity = 1;

	template<typename NumType>
	inline
	MaxTimes<NumType>::MaxTimes(const NumType value)
		: value(value) {
	}

	template<typename NumType>
	inline MaxTimes<NumType>&
	MaxTimes<NumType>::operator*=(const MaxTimes<NumType>& x) {
		value *= x.value;
		return *this;
	}

	template<typename NumType>
	inline MaxTimes<NumType>&
	MaxTimes<NumType>::operator+=(const MaxTimes<NumType>& x) {
		value = std::max(value, x.value);
		return *this;
	}

	template<typename NumType>
	inline MaxTimes<NumType>
	MaxTimes<NumType>::operator*(const MaxTimes<NumType>& x) const {
		return MaxTimes<NumType>(value * x.value);
	}

	template<typename NumType>
	inline MaxTimes<NumType>
	MaxTimes<NumType>::operator+(const MaxTimes<NumType>& x) const {
		return MaxTimes<NumType>(std::max(value, x.value));
	}

	template<typename NumType>
	inline bool
	MaxTimes<NumType>::operator<(const MaxTimes<NumType>& x) const {
		return value < x.value;
	}

	template<typename NumType>
	inline bool
	MaxTimes<NumType>::operator==(const MaxTimes<NumType>& x) const {
		return value == x.value;
	}
		
}
		
#endif // __MAX_TIMES_HH__
