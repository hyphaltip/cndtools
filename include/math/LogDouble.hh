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

#ifndef __MATH_LOGDOUBLE_HH__
#define __MATH_LOGDOUBLE_HH__

#include <stdexcept>
#include <cmath>
#include <istream>
#include <ostream>

#include "PiecewiseLinearInterpolator.hh"
#include "PiecewiseCubicHermiteInterpolator.hh"

namespace math {

	class LogDouble {
	public:
		LogDouble(double val = 0.0);

		LogDouble& operator*=(const LogDouble& other);
		LogDouble& operator+=(const LogDouble& other);
		LogDouble& operator/=(const LogDouble& other);
		LogDouble& operator-=(const LogDouble& other);

		LogDouble operator*(const LogDouble& other) const;
		LogDouble operator+(const LogDouble& other) const;
		LogDouble operator/(const LogDouble& other) const;
		LogDouble operator-(const LogDouble& other) const;

		LogDouble& operator^=(const double& exponent);
		LogDouble operator^(const double& exponent) const;

		bool operator==(const LogDouble& other) const;
		bool operator>(const LogDouble& other) const;
		bool operator<(const LogDouble& other) const;
		bool operator<=(const LogDouble& other) const;
		bool operator>=(const LogDouble& other) const;

		operator double() const { return std::exp(value); }

		double value;

	private:
		static double log_one_plus_exp(double x);
		static double log_one_minus_exp(double x);
		static double derivative_log_one_plus_exp(double x);
		static PiecewiseLinearInterpolator log_one_plus_exp_approx;
		static PiecewiseCubicHermiteInterpolator log_one_plus_exp_approx2;
	};

	inline LogDouble::LogDouble(double val)
		: value(std::log(val)) {
	}

	inline LogDouble& LogDouble::operator*=(const LogDouble& other) {
		value += other.value;
		return *this;
	}

	inline LogDouble& LogDouble::operator/=(const LogDouble& other) {
		value -= other.value;
		return *this;
	}

	inline LogDouble& LogDouble::operator^=(const double& exponent) {
		value *= exponent;
		return *this;
	}
	
	inline LogDouble LogDouble::operator*(const LogDouble& other) const {
		return LogDouble(*this) *= other;
	}

	inline LogDouble LogDouble::operator+(const LogDouble& other) const {
		return LogDouble(*this) += other;
	}

	inline LogDouble LogDouble::operator/(const LogDouble& other) const {
		return LogDouble(*this) /= other;
	}

	inline LogDouble LogDouble::operator-(const LogDouble& other) const {
		return LogDouble(*this) -= other;
	}	

	inline LogDouble LogDouble::operator^(const double& exponent) const {
		return LogDouble(*this) ^= exponent;
	}
	
	inline bool LogDouble::operator==(const LogDouble& other) const {
		return value == other.value;
	}

	inline bool LogDouble::operator<(const LogDouble& other) const {
		return value < other.value;
	}

	inline bool LogDouble::operator>(const LogDouble& other) const {
		return value > other.value;
	}

	inline bool LogDouble::operator>=(const LogDouble& other) const {
		return value >= other.value;
	}

	inline bool LogDouble::operator<=(const LogDouble& other) const {
		return value <= other.value;
	}

	inline std::istream& operator>>(std::istream& stream, LogDouble& ld) {
		return stream >> ld.value;
	}

	inline std::ostream& operator<<(std::ostream& stream, const LogDouble& ld) {
		return stream << ld.value;
	}
	
}
		
#endif // __MATH_LOGDOUBLE_HH__
