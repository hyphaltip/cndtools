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

#ifndef __PARSER_HH__
#define __PARSER_HH__

#include <string>
#include <istream>
#include <limits>
#include "boost/tokenizer.hpp"

namespace util {

	namespace parser {

		// Exception raised when an invalid format is found
		class FormatError {
		private:
			std::string line;
			std::string problem;
		public:
			FormatError(const std::string& line, const std::string& problem) :
				line(line),
				problem(problem)
			{}
			const std::string& getLine() const throw() { return line; }
			const std::string& getProblem() const throw() { return problem; }
		};

		const std::string numericChars = "1234567890.";
	
		inline bool isFreeText(const std::string& s) {
			return s.find_first_not_of(numericChars) != std::string::npos;
		}

		inline void ignoreLine(std::istream& strm) {
			strm.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		inline void ignoreWhitespace(std::istream& strm) {
			while (!strm.eof() && std::isspace(strm.peek())) {
				strm.ignore();
			}
		}

		// Split the string S at each tab and return the tokens in vector V
		inline void splitOnTabs(const std::string& s,
								std::vector<std::string>& v) {
			typedef boost::tokenizer<boost::char_separator<char> > CharTok;
			CharTok tabtok(s, boost::char_separator<char>("\t"));
			std::copy(tabtok.begin(), tabtok.end(), std::back_inserter(v));
		}

	};

};

#endif // __PARSER_HH__
