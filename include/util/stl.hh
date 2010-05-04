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

#ifndef __STL_HH__
#define __STL_HH__

#include <functional>
#include <numeric>
#include <ostream>
#include <algorithm>
#include <string>

namespace util {

	namespace stl {
    
		template<typename InputIterator1, typename InputIterator2>
		size_t matches(InputIterator1 first1, InputIterator1 last1,
					   InputIterator2 first2) {
			typedef typename std::iterator_traits<InputIterator1>::value_type
				value_type;
			return std::inner_product(first1, last1, first2, 0,
									  std::plus<size_t>(),
									  std::equal_to<value_type>());
		}

		template<typename InputIterator1, typename InputIterator2>
		size_t mismatches(InputIterator1 first1, InputIterator1 last1,
						  InputIterator2 first2) {
			typedef typename std::iterator_traits<InputIterator1>::value_type
				value_type;
			return std::inner_product(first1, last1, first2, 0,
									  std::plus<size_t>(),
									  std::not_equal_to<value_type>());
		}

		template<typename Iterator>
		void print_elements(std::ostream& strm,
						   Iterator begin,
						   Iterator end,
						   const std::string& sep = " ") {
			Iterator pos = begin;
			while (pos != end) {
				if (pos != begin) {
					strm << sep;
				}
				strm << *pos;
				++pos;
			}
			strm << '\n';
		}

		template<typename Container>
		void print_elements(std::ostream& strm,
							const Container& c,
							const std::string& sep = " ") {
			print_elements(strm, c.begin(), c.end(), sep);
		}
		
		template<typename InputIterator>
		typename InputIterator::value_type
		sum(InputIterator begin,
			InputIterator end,
			typename InputIterator::value_type init =
			typename InputIterator::value_type())
		{
			while (begin != end) {
				init += *begin;
				++begin;
			}
			return init;
		}

		template<typename InputIterator>
		typename InputIterator::value_type
		product(InputIterator begin,
				InputIterator end,
				typename InputIterator::value_type init =
				typename InputIterator::value_type())
		{
			while (begin != end) {
				init *= *begin;
				++begin;
			}
			return init;
		}
		
		template<typename InputIterator, typename UnaryFunc>
		typename UnaryFunc::result_type
		sum(InputIterator begin,
			InputIterator end,
			UnaryFunc f,
			typename UnaryFunc::result_type initValue = typename UnaryFunc::result_type())
		{
			while (begin != end) {
				initValue += f(*begin);
				++begin;
			}
			return initValue;
		}

		template<typename InputIterator, typename UnaryFunc>
		typename UnaryFunc::result_type
		product(InputIterator begin,
				InputIterator end,
				UnaryFunc f,
				typename UnaryFunc::result_type initValue = typename UnaryFunc::result_type())
		{
			while (begin != end) {
				initValue *= f(*begin);
				++begin;
			}
			return initValue;
		}
		
		template<class OP>
		struct FunctorRef : public std::unary_function<typename OP::argument_type,
													   typename OP::result_type> {
			const OP& op;
			FunctorRef(const OP& op) : op(op) {}
			typename OP::result_type
			operator()(typename OP::argument_type x) {
				return op(x);
			}
		};

		template<class OP>
		FunctorRef<OP> make_functor_ref(const OP& op) { return FunctorRef<OP>(op); }

		template<typename Container>
		void removeDuplicates(Container& c) {
			std::sort(c.begin(), c.end());
			c.erase(std::unique(c.begin(), c.end()), c.end());
		}

		// Implementation of copy_if from Effective STL
		template<typename InputIterator,
				 typename OutputIterator,
				 typename Predicate>
		OutputIterator copy_if(InputIterator begin,
							   InputIterator end,
							   OutputIterator destBegin,
							   Predicate p)
		{
			while (begin != end) {
				if (p(*begin)) *destBegin++ = *begin;
				++begin;
			}
			return destBegin;
		}

		// Generic iterator over an input stream
		template<typename StreamType>
		class InputStreamIterator
			: public std::iterator<std::input_iterator_tag,
								   typename StreamType::ValueType>
		{
		public:
			typedef typename StreamType::ValueType ValueType;		

		private:
			StreamType* strm;
			ValueType value;
			bool ok;
		
		public:      
			InputStreamIterator() : strm(NULL), ok(false) {
			}

			InputStreamIterator(StreamType& s) : strm(&s) {
				read();
			}

			InputStreamIterator(const InputStreamIterator& obj) 
				: strm(obj.strm), value(obj.value), ok(obj.ok) {
			}
		
			ValueType& operator*() {
				return value;
			}
		
			ValueType* operator->() {
				return &(operator*());
			}

			InputStreamIterator& operator++() {
				read(); 
				return *this; 
			}

			InputStreamIterator operator++(int)  {
				InputStreamIterator tmp = *this;
				read();
				return tmp;
			}

			bool equals(const InputStreamIterator& x) const {
				return (ok == x.ok) && (!ok || strm == x.strm);
			}
		
		private:      

			void read() {
				ok = (strm && *strm) ? true : false;
				if (ok) {
					*strm >> value;
					ok = *strm ? true : false;
				}
			}
		};
  
		template<typename StreamType>
		inline bool 
		operator==(const InputStreamIterator<StreamType>& x,
				   const InputStreamIterator<StreamType>& y) {
			return x.equals(y);
		}

		template<typename StreamType>
		inline bool 
		operator!=(const InputStreamIterator<StreamType>& x,
				   const InputStreamIterator<StreamType>& y) {
			return !x.equals(y);
		}
	
	};
	
};

#endif // __STL_HH__
