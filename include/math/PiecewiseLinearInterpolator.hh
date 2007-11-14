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

#ifndef __MATH_PIECEWISELINEARINTERPOLATOR_HH__
#define __MATH_PIECEWISELINEARINTERPOLATOR_HH__

#include <vector>
#include <cmath>

namespace math {

    class PiecewiseLinearInterpolator {
	public:
		template<typename Function>
		PiecewiseLinearInterpolator(Function f,
									double min,
									double max,
									int log_sample_rate);

		double operator()(double x);
		
	private:
		const double min;
		const double max;
		const int log_sample_rate;
		const double sample_rate;
		std::vector<double> y;
    };

	template<typename Function>
	PiecewiseLinearInterpolator::
	PiecewiseLinearInterpolator(Function f,
								double min,
								double max,
								int log_sample_rate) :
		min(min),
		max(max),
		log_sample_rate(log_sample_rate),
		sample_rate(std::ldexp(1.0, log_sample_rate)),
		y(static_cast<size_t>(std::ldexp(max - min, log_sample_rate)) + 2, 0.0)
	{
		for (size_t i = 0; i < y.size(); ++i) {
			double x = min + std::ldexp(static_cast<double>(i),
										-log_sample_rate);
			y[i] = f(x);
		}
	}

	inline double PiecewiseLinearInterpolator::operator()(double x) {
		if (x < min or x > max) {
			return 0.0;
		} else {
			double i = (x - min) * sample_rate;
			int k = static_cast<int>(i);
			double s = i - k;
			return y[k] + s * (y[k + 1] - y[k]);
		}
	}
	
}

#endif // __MATH_PIECEWISELINEARINTERPOLATOR_HH__
