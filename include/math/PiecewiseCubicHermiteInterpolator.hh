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

#ifndef __MATH_PIECEWISECUBICHERMITEINTERPOLATOR_HH__
#define __MATH_PIECEWISECUBICHERMITEINTERPOLATOR_HH__

#include <vector>
#include <cmath>

namespace math {

	class PiecewiseCubicHermiteInterpolator {
	public:
		template<typename Function, typename Derivative>
		PiecewiseCubicHermiteInterpolator(Function f,
										  Derivative df,
										  double min,
										  double max,
										  int log_sample_rate);
		
		double operator()(double x);

	private:
		struct HermiteValues {
			double y;
			double d;
			double c;
			double b;
		};
			
		const double min;
		const double max;
		const int log_sample_rate;
		const double sample_rate;
		const double sample_rate_inverse;
		std::vector<HermiteValues> values;
	};
	
	template<typename Function, typename Derivative>
	PiecewiseCubicHermiteInterpolator::
	PiecewiseCubicHermiteInterpolator(Function f,
									  Derivative df,
									  double min,
									  double max,
									  int log_sample_rate) :
		min(min),
		max(max),
		log_sample_rate(log_sample_rate),
		sample_rate(std::ldexp(1.0, log_sample_rate)),
		sample_rate_inverse(std::ldexp(1.0, -log_sample_rate)),
		values(static_cast<size_t>(std::ldexp(max - min, log_sample_rate)) + 2)
	{
		for (size_t i = 0; i < values.size(); ++i) {
				double x = min + std::ldexp(static_cast<double>(i),
											-log_sample_rate);
				values[i].y = f(x);
				values[i].d = df(x);
			}
				 for (size_t i = 0; i < (values.size() - 1); ++i) {
					 double delta = std::ldexp(values[i + 1].y - values[i].y, log_sample_rate);
					 values[i].c = std::ldexp(3 * delta - 2 * values[i].d - values[i + 1].d, log_sample_rate);
					 values[i].b = std::ldexp(values[i].d - 2 * delta + values[i + 1].d, log_sample_rate * 2);
				 }
				 }
		
		inline double PiecewiseCubicHermiteInterpolator::operator()(double x) {
			if (x < min or x > max) {
				return 0.0;
			} else {
				double i = (x - min) * sample_rate;
				int k = static_cast<int>(i);
				HermiteValues v(values[k]);
				double s = (i - k) * sample_rate_inverse;
				return v.y + s * (v.d + s * (v.c + s * v.b));
			}
		}
	
	}

#endif // __MATH_PIECEWISECUBICHERMITEINTERPOLATOR_HH__
