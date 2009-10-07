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

#ifndef __UTIL_MATRIX_HH__
#define __UTIL_MATRIX_HH__

#include <valarray>
#include <ostream>

namespace util {

	template<typename NumType>
	class Matrix {
	public:
		Matrix(size_t nRows = 0, size_t nCols = 0);
		Matrix(size_t nRows, size_t nCols, NumType val);	

		size_t size() const;
		size_t getNumRows() const;
		size_t getNumCols() const;

		void transpose();
		void reverse();
		void reverseRows();
		void reverseCols();
		
		void resize(size_t nRows, size_t nCols, NumType val = NumType());
		void fill(NumType val);
	
		NumType& operator()(size_t i, size_t j);
		NumType operator()(size_t i, size_t j) const;

		Matrix& operator+=(const Matrix& m);
		Matrix& operator-=(const Matrix& m);
		Matrix& operator*=(const Matrix& m);
		Matrix& operator/=(const Matrix& m);
		Matrix& operator%=(const Matrix& m);
		Matrix& operator^=(const Matrix& m);
		Matrix& operator&=(const Matrix& m);
		Matrix& operator|=(const Matrix& m);
		Matrix& operator<<=(const Matrix& m);
		Matrix& operator>>=(const Matrix& m);

		Matrix& operator+=(const NumType& x);
		Matrix& operator-=(const NumType& x);
		Matrix& operator*=(const NumType& x);
		Matrix& operator/=(const NumType& x);
		Matrix& operator%=(const NumType& x);
		Matrix& operator^=(const NumType& x);
		Matrix& operator&=(const NumType& x);
		Matrix& operator|=(const NumType& x);
		Matrix& operator<<=(const NumType& x);
		Matrix& operator>>=(const NumType& x);

	private:
		size_t nRows;
		size_t nCols;
		std::valarray<NumType> v;
	};

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator+=(const Matrix& m) {
		v += m.v;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator-=(const Matrix& m) {
		v -= m.v;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator*=(const Matrix& m) {
		v *= m.v;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator/=(const Matrix& m) {
		v /= m.v;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator%=(const Matrix& m) {
		v %= m.v;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator^=(const Matrix& m) {
		v ^= m.v;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator&=(const Matrix& m) {
		v &= m.v;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator|=(const Matrix& m) {
		v |= m.v;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator<<=(const Matrix& m) {
		v <<= m.v;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator>>=(const Matrix& m) {
		v >>= m.v;
		return *this;
	}
	

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator+=(const NumType& x) {
		v += x;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator-=(const NumType& x) {
		v -= x;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator*=(const NumType& x) {
		v *= x;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator/=(const NumType& x) {
		v /= x;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator%=(const NumType& x) {
		v %= x;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator^=(const NumType& x) {
		v ^= x;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator&=(const NumType& x) {
		v &= x;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator|=(const NumType& x) {
		v |= x;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator<<=(const NumType& x) {
		v <<= x;
		return *this;
	}

	template<typename NumType>
	inline Matrix<NumType>& Matrix<NumType>::operator>>=(const NumType& x) {
		v >>= x;
		return *this;
	}

	template<typename NumType>
	void Matrix<NumType>::transpose() {
		Matrix<NumType> old(*this);
		resize(getNumCols(), getNumRows());
		for (size_t i = 0; i < getNumRows(); ++i) {
			for (size_t j = 0; j < getNumCols(); ++j) {
				(*this)(i, j) = old(j, i);
			}
		}
	}

	template<typename NumType>
	void Matrix<NumType>::reverse() {
		Matrix<NumType> old(*this);
		for (size_t i = 0; i < getNumRows(); ++i) {
			for (size_t j = 0; j < getNumCols(); ++j) {
				(*this)(i, j) = old(getNumRows() - i - 1,
									getNumCols() - j - 1);
			}
		}
	}

	template<typename NumType>
	void Matrix<NumType>::reverseRows() {
		Matrix<NumType> old(*this);
		for (size_t i = 0; i < getNumRows(); ++i) {
			for (size_t j = 0; j < getNumCols(); ++j) {
				(*this)(i, j) = old(getNumRows() - i - 1, j);
			}
		}
	}
	
	template<typename NumType>
	void Matrix<NumType>::reverseCols() {
		Matrix<NumType> old(*this);
		for (size_t i = 0; i < getNumRows(); ++i) {
			for (size_t j = 0; j < getNumCols(); ++j) {
				(*this)(i, j) = old(i, getNumCols() - j - 1);
			}
		}
	}

	template<typename NumType>
	std::ostream& operator<<(std::ostream& strm, const Matrix<NumType>& m) {
		for (size_t i = 0; i < m.getNumRows(); ++i) {
			for (size_t j = 0; j < m.getNumCols(); ++j) {
				if (j != 0) {
					strm << '\t';
				}
				strm << m(i, j);
			}
			strm << '\n';
		}
		return strm;
	}

	template<typename NumType>
	Matrix<NumType>::Matrix(size_t nRows, size_t nCols)
		: nRows(nRows), nCols(nCols), v(nRows * nCols) {
	}

	template<typename NumType>
	Matrix<NumType>::Matrix(size_t nRows, size_t nCols, NumType val)
		: nRows(nRows), nCols(nCols), v(val, nRows * nCols) {
	}

	template<typename NumType>
	void Matrix<NumType>::resize(size_t nRows, size_t nCols, NumType val) {
		size_t oldSize = this->nRows * this->nCols;
		this->nRows = nRows;
		this->nCols = nCols;
		if (oldSize != nRows * nCols) {
			v.resize(nRows * nCols, val);
		}
	}

	template<typename NumType>
	void Matrix<NumType>::fill(NumType val) {
		for (size_t i = 0; i < getNumRows(); ++i) {
			for (size_t j = 0; j < getNumCols(); ++j) {
				(*this)(i, j) = val;
			}
		}
	}
	
	template<typename NumType>
	inline size_t Matrix<NumType>::size() const { return v.size(); }

	template<typename NumType>
	inline size_t Matrix<NumType>::getNumRows() const { return nRows; }

	template<typename NumType>
	inline size_t Matrix<NumType>::getNumCols() const { return nCols; }

	template<typename NumType>
	inline NumType& Matrix<NumType>::operator()(size_t i, size_t j) {
		return v[i * nCols + j];
	}

	template<typename NumType>
	inline NumType Matrix<NumType>::operator()(size_t i, size_t j) const {
		return v[i * nCols + j];
	}

};
	
#endif // __UTIL_MATRIX_HH__
