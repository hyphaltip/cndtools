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

#ifndef __PLUS_TIMES_HH__
#define __PLUS_TIMES_HH__

namespace math {

	template<typename NumType>
	class PlusTimes {
	public:
		typedef PlusTimes<NumType> Element;
		static const NumType zero;
		static const NumType multiplicativeIdentity;

		NumType value;

		PlusTimes<NumType>(const NumType value = NumType());
	
		PlusTimes<NumType>& operator*=(const PlusTimes<NumType>& x);
		PlusTimes<NumType>& operator+=(const PlusTimes<NumType>& x);

		PlusTimes<NumType> operator*(const PlusTimes<NumType>& x) const;
		PlusTimes<NumType> operator+(const PlusTimes<NumType>& x) const;

		bool operator<(const PlusTimes<NumType>& x) const;	
		bool operator==(const PlusTimes<NumType>& x) const;
	};
		
	template<typename NumType>
	const NumType PlusTimes<NumType>::zero = 0;

	template<typename NumType>
	const NumType PlusTimes<NumType>::multiplicativeIdentity = 1;

	template<typename NumType>
	inline
	PlusTimes<NumType>::PlusTimes(const NumType value)
		: value(value) {
	}

	template<typename NumType>
	inline PlusTimes<NumType>&
	PlusTimes<NumType>::operator*=(const PlusTimes<NumType>& x) {
		value *= x.value;
		return *this;
	}

	template<typename NumType>
	inline PlusTimes<NumType>&
	PlusTimes<NumType>::operator+=(const PlusTimes<NumType>& x) {
		value += x.value;
		return *this;
	}

	template<typename NumType>
	inline PlusTimes<NumType>
	PlusTimes<NumType>::operator*(const PlusTimes<NumType>& x) const {
		return PlusTimes<NumType>(value * x.value);
	}

	template<typename NumType>
	inline PlusTimes<NumType>
	PlusTimes<NumType>::operator+(const PlusTimes<NumType>& x) const {
		return PlusTimes<NumType>(value + x.value);
	}

	template<typename NumType>
	inline bool
	PlusTimes<NumType>::operator<(const PlusTimes<NumType>& x) const {
		return value < x.value;
	}

	template<typename NumType>
	inline bool
	PlusTimes<NumType>::operator==(const PlusTimes<NumType>& x) const {
		return value == x.value;
	}
		
}
		
#endif // __PLUS_TIMES_HH__
