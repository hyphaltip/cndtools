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

#include <limits>

#include "math/LogDouble.hh"

namespace math {

	double LogDouble::log_one_plus_exp(double x) {
		return std::log(1.0 + std::exp(x));
	}

	double LogDouble::log_one_minus_exp(double x) {
		return std::log(1.0 - std::exp(x));
	}
	
	double LogDouble::derivative_log_one_plus_exp(double x) {
		return std::exp(x) / (1.0 + std::exp(x));
	}
	
	PiecewiseLinearInterpolator
	LogDouble::log_one_plus_exp_approx(log_one_plus_exp, -37.0, 0.0, 10);

	PiecewiseCubicHermiteInterpolator
	LogDouble::log_one_plus_exp_approx2(log_one_plus_exp, derivative_log_one_plus_exp, -37.0, 0.0, 8);
	
	LogDouble& LogDouble::operator+=(const LogDouble& other) {
		double other_value(other.value);
		if (value < other_value) { std::swap(value, other_value); }
		if (other_value != -std::numeric_limits<double>::infinity()) { 
			value += log_one_plus_exp_approx(other_value - value);
		}
		return *this;
	}
	
	LogDouble& LogDouble::operator-=(const LogDouble& other) {
		if (value < other.value) {
			throw std::runtime_error("Attempted to subtract large value "
									 "from small value in log space");
		} else if (other.value != -std::numeric_limits<double>::infinity()) {
			value += log_one_minus_exp(other.value - value);
		}
		return *this;
	}

}
