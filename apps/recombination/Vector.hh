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

#ifndef __VECTOR_HH__
#define __VECTOR_HH__

#include <valarray>
#include <iosfwd>

template<typename T>
class Vector {
public:
	static Vector<T> unit(size_t dims, size_t dim);
		
	Vector(size_t dims = 0);

	Vector& operator=(const Vector& x);
		
	Vector& operator+=(const Vector& x);
	Vector& operator+=(const T x);
	Vector& operator-=(const Vector& x);
	Vector& operator-=(const T x);
	Vector& operator*=(const Vector& x);
	Vector& operator*=(const T x);
	Vector& operator/=(const Vector& x);
	Vector& operator/=(const T x);

	Vector operator+(const Vector& x) const;
	Vector operator+(const T x) const;
	Vector operator-(const Vector& x) const;
	Vector operator-(const T x) const;
	Vector operator*(const Vector& x) const;
	Vector operator*(const T x) const;
	Vector operator/(const Vector& x) const;
	Vector operator/(const T x) const;

	bool operator<(const Vector& x) const;
	bool operator==(const Vector& x) const;

	T& operator[](size_t i);
	const T& operator[](size_t i) const;

	size_t size() const;

	template<typename S>
	friend std::ostream& operator<<(std::ostream& stream, const Vector<S>& v);

	template<typename S>
	friend std::istream& operator>>(std::istream& stream, Vector<S>& v);
	
protected:
	std::valarray<T> v;
};

template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& x) {
	if (this == &x) { return *this; }
	if (v.size() != x.v.size()) { v.resize(x.v.size()); }
	v = x.v;
	return *this;
}
	
template<typename T>
Vector<T> Vector<T>::unit(size_t dims, size_t dim) {
	Vector<T> v(dims);
	v[dim] = 1;
	return v;
}
	
template<typename T>
Vector<T>::Vector(size_t dims) : v(dims) {}

template<typename T>
Vector<T>& Vector<T>::operator+=(const Vector& x) {
	v += x.v;
	return *this;
}
	
template<typename T>
Vector<T>& Vector<T>::operator+=(const T x) {
	v += x;
	return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator-=(const Vector& x) {
	v -= x.v;
	return *this;
}
	
template<typename T>
Vector<T>& Vector<T>::operator-=(const T x) {
	v -= x;
	return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator*=(const Vector& x) {
	v *= x.v;
	return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator*=(const T x) {
	v *= x;
	return *this;
}
	
template<typename T>
Vector<T>& Vector<T>::operator/=(const Vector& x) {
	v /= x.v;
	return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator/=(const T x) {
	v /= x;
	return *this;
}

template<typename T>
Vector<T> Vector<T>::operator+(const Vector& x) const {
	return Vector<T>(*this) += x;
}
	
template<typename T>
Vector<T> Vector<T>::operator+(const T x) const {
	return Vector<T>(*this) += x;
}

template<typename T>
Vector<T> Vector<T>::operator-(const Vector& x) const {
	return Vector<T>(*this) -= x;
}
	
template<typename T>
Vector<T> Vector<T>::operator-(const T x) const {
	return Vector<T>(*this) -= x;
}

template<typename T>
Vector<T> Vector<T>::operator*(const Vector& x) const {
	return Vector<T>(*this) *= x;
}
	
template<typename T>
Vector<T> Vector<T>::operator*(const T x) const {
	return Vector<T>(*this) *= x;
}

template<typename T>
Vector<T> Vector<T>::operator/(const Vector& x) const {
	return Vector<T>(*this) /= x;
}
	
template<typename T>
Vector<T> Vector<T>::operator/(const T x) const {
	return Vector<T>(*this) /= x;
}

template<typename T>
bool
Vector<T>::operator<(const Vector& x) const {
	for (size_t i = 0; i < v.size(); ++i) {
		if (v[i] < x.v[i]) {
			return true;
		} else if (x.v[i] < v[i]) {
			return false;
		}
	}
	return false;
}

template<typename T>
bool
Vector<T>::operator==(const Vector& x) const {
	for (size_t i = 0; i < v.size(); ++i) {
		if (v[i] != x.v[i]) {
			return false;
		}
	}
	return true;
}
	
template<typename T>
inline T& Vector<T>::operator[](size_t i) { return v[i]; }
	
template<typename T>
inline const T& Vector<T>::operator[](size_t i) const { return v[i]; }

template<typename T>
inline size_t Vector<T>::size() const { return v.size(); }

template<typename T>
std::ostream& operator<<(std::ostream& stream, const Vector<T>& v) {
	for (size_t i = 0; i < v.v.size(); ++i) {
		if (i > 0) { stream << '\t'; }
		stream << v.v[i];
	}
	return stream;
}

template<typename T>
std::istream& operator>>(std::istream& stream, Vector<T>& v) {
	for (size_t i = 0; i < v.v.size(); ++i) {
		stream >> v.v[i];
	}
	return stream;
}

#endif // __VECTOR_HH__
