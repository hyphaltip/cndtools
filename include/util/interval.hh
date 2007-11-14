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

#ifndef __INTERVAL_HH__
#define __INTERVAL_HH__

#include <queue>
#include <algorithm>
#include <functional>
#include <utility>

namespace util {

	namespace interval {

		template<class T>
		struct interval_traits {
			typedef typename T::coord_type coord_type;
			typedef typename T::difference_type difference_type;
			static coord_type start(const T& i) { return i.start(); }
			static coord_type end(const T& i) { return i.end(); }
		};

		template<class T>
		struct interval_traits<T*> {
			typedef typename T::coord_type coord_type;
			typedef typename T::difference_type difference_type;
			static coord_type start(const T* i) { return i->start(); }
			static coord_type end(const T* i) { return i->end(); }
		};
	
		template<class T1, class T2>
		bool intervalsOverlap(const T1& i1, const T2& i2) {
			typedef interval_traits<T1> IT1;
			typedef interval_traits<T2> IT2;
			return (IT1::start(i1) < IT2::end(i2) &&
					IT1::end(i1) > IT2::start(i2));
		}
	
		template<class T1, class T2>
		typename interval_traits<T1>::difference_type
		overlapLength(const T1& i1, const T2& i2) {
			typedef interval_traits<T1> IT1;
			typedef interval_traits<T2> IT2;
			return std::min(IT1::end(i1), IT2::end(i2))
				- std::max(IT1::start(i1), IT2::start(i2));
		}
	
		template<class T1, class T2>
		typename interval_traits<T1>::difference_type
		distanceBetween(const T1& i1, const T2& i2) {
			typedef interval_traits<T1> IT1;
			typedef interval_traits<T2> IT2;		
			return std::max(IT1::start(i1), IT2::start(i2))
				- std::min(IT1::end(i1), IT2::end(i2));
		}
	
		template<class T>
		struct StartSorter {
			typedef interval_traits<T> IT;
			bool operator()(const T& i1, const T& i2) {
				return IT::start(i1) < IT::start(i2);
			}
		};

		template<class T, class IT = interval_traits<T> >
		struct SmallestEndFirst {
			bool operator()(const T& i1, const T& i2) {
				return IT::end(i2) < IT::end(i1);
			}
		};

		template<typename InputIterator, typename OutputIterator>
		void overlaps(InputIterator first,
					  InputIterator last,
					  OutputIterator dest) {
			typedef typename InputIterator::value_type T;
			typedef interval_traits<T> IT;
		
			// The intervals that we are currently inside
			std::vector<T> inside;
	
			// Loop until we have passed through all intervals
			while (first != last) {
				// Check if we currently inside any intervals
				if (inside.empty()) {
					inside.push_back(*first);
					++first;
				}
				// Check if the next interval does not overlap with the
				// first of the intervals that we are currently inside
				else if (IT::end(inside.front()) <= IT::start(*first)) {
					std::pop_heap(inside.begin(), inside.end(),
								  SmallestEndFirst<T>());
					inside.pop_back();
				}
				// Otherwise, the next interval overlaps with all of the
				// intervals currently in inside
				else {
					// Output pairs of the next interval with all of the
					// intervals we are inside
					std::transform(inside.begin(), inside.end(), dest,
								   std::bind2nd(std::ptr_fun(&std::make_pair<T, T>),
												*first));
					inside.push_back(*first);
					std::make_heap(inside.begin(), inside.end(),
								   SmallestEndFirst<T>());
					++first;
				}
			}
		}


		template<typename InputIterator1,
				 typename InputIterator2,
				 typename OutputIterator>
		void overlaps(InputIterator1 first1, InputIterator1 last1,
					  InputIterator2 first2, InputIterator2 last2,
					  OutputIterator dest) {
			interval_traits<typename InputIterator1::value_type> it1;
			interval_traits<typename InputIterator2::value_type> it2;
			return overlaps(first1, last1, first2, last2, dest, it1, it2);
		}

		template<typename InputIterator1,
				 typename InputIterator2,
				 typename OutputIterator,
				 typename IntervalTraits1,
				 typename IntervalTraits2>
		void overlaps(InputIterator1 first1, InputIterator1 last1,
					  InputIterator2 first2, InputIterator2 last2,
					  OutputIterator dest,
					  const IntervalTraits1& it1,
					  const IntervalTraits2& it2) {
			typedef typename InputIterator1::value_type T1;
			typedef typename InputIterator2::value_type T2;
			
			// The interval that we are currently inside
			std::vector<T1> inside1;
			std::vector<T2> inside2;
			
			// Loop until we have passed through all intervals
			while (first1 != last1 or first2 != last2) {
				if ((first2 == last2) or
					((first1 != last1) and
					 it1.start(*first1) < it2.start(*first2))) {
					// Check if the next interval does not overlap with
					// the first of the intervals that we are currently
					// inside for the first sequence
					if (!inside1.empty() &&
						it1.end(inside1.front()) <= it1.start(*first1)) {
						std::pop_heap(inside1.begin(), inside1.end(),
									  SmallestEndFirst<T1, IntervalTraits1>());
						inside1.pop_back();
						// Check if the next interval does not overlap with
						// the first of the intervals that we are currently
						// inside for the second sequence
					} else if (!inside2.empty() &&
							   it2.end(inside2.front()) <= it1.start(*first1)) {
						std::pop_heap(inside2.begin(), inside2.end(),
									  SmallestEndFirst<T2, IntervalTraits2>());
						inside2.pop_back();
					} else {
						// Output pairs of the next interval with all of the
						// intervals we are inside
						std::transform(inside2.begin(), inside2.end(), dest,
									   std::bind1st(std::ptr_fun(&std::make_pair<T1, T2>),
													*first1));
						inside1.push_back(*first1);
						std::make_heap(inside1.begin(), inside1.end(),
									   SmallestEndFirst<T1, IntervalTraits1>());
						++first1;
					}
				} else {
					// Check if the next interval does not overlap with
					// the first of the intervals that we are currently
					// inside for the first sequence
					if (!inside1.empty() &&
						it1.end(inside1.front()) <= it2.start(*first2)) {
						std::pop_heap(inside1.begin(), inside1.end(),
									  SmallestEndFirst<T1, IntervalTraits1>());
						inside1.pop_back();
						// Check if the next interval does not overlap with
						// the first of the intervals that we are currently
						// inside for the second sequence
					} else if (!inside2.empty() &&
							   it2.end(inside2.front()) <= it2.start(*first2)) {
						std::pop_heap(inside2.begin(), inside2.end(),
									  SmallestEndFirst<T2, IntervalTraits2>());
						inside2.pop_back();
					} else {
						// Output pairs of the next interval with all of the
						// intervals we are inside
						std::transform(inside1.begin(), inside1.end(), dest,
									   std::bind2nd(std::ptr_fun(&std::make_pair<T1, T2>),
													*first2));
						inside2.push_back(*first2);
						std::make_heap(inside2.begin(), inside2.end(),
									   SmallestEndFirst<T2, IntervalTraits2>());
						++first2;
					}
				}
			}
		}

	};

	
};

#endif // __INTERVAL_HH__
