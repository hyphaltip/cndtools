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

#ifndef __CHARMANIP_HH__
#define __CHARMANIP_HH__

#include <string>
#include <cctype>

namespace util {

	namespace charmanip {

		struct Mapper : public std::unary_function<const char, char> {
			char table[256];
			Mapper(const std::string from,
				   const std::string to,
				   bool caseSensitive=true) {
				for (unsigned int i = 0; i < 256; ++i) {
					table[i] = i;
				}
				for (unsigned int i = 0; i < from.size(); ++i) {
					if (caseSensitive) {
						table[static_cast<unsigned char>(from[i])] = to[i];
					} else {
						table[static_cast<unsigned char>(std::toupper(from[i]))] = std::toupper(to[i]);
						table[static_cast<unsigned char>(std::tolower(from[i]))] = std::tolower(to[i]);
					}					
				}
			}
			char operator()(const char c) const {
				return table[static_cast<unsigned char>(c)];
			}
		};
	
		class Encoder : public std::unary_function<const char, unsigned char> {
		private:
			unsigned char table[256];		

			void init(const std::string& chars, bool caseSensitive) {
				for (unsigned int i = 0; i < 256; ++i) {
					table[i] = 0xFF;
				}
				for (unsigned int i = 0; i < chars.size(); ++i) {
					if (caseSensitive) {
						table[static_cast<unsigned char>(chars[i])] = i;
					} else {
						table[static_cast<unsigned char>(std::toupper(chars[i]))] = i;
						table[static_cast<unsigned char>(std::tolower(chars[i]))] = i;
					}
				}
			}
		public:

			Encoder(const std::string chars,
					bool caseSensitive=true) {
				init(chars, caseSensitive);
			}

			Encoder(const std::string chars,
					bool caseSensitive,
					char nonMemberChar) {
			
				init(chars, caseSensitive);
				unsigned char nonMemberUpperEncoding =
					table[static_cast<unsigned int>(std::toupper(nonMemberChar))];
				unsigned char nonMemberLowerEncoding =
					table[static_cast<unsigned int>(std::tolower(nonMemberChar))];
				for (unsigned int i = 0; i < 256; ++i) {
					if (table[i] == 0xFF) {
						if (std::isupper(static_cast<char>(i))) {
							table[i] = nonMemberUpperEncoding;
						} else {
							table[i] = nonMemberLowerEncoding;
						}
					}
				}
			}
		
			unsigned char operator()(const char c) const {
				return table[static_cast<unsigned char>(c)];
			}
			bool isMember(const char c) const {
				return table[static_cast<unsigned char>(c)] != 0xFF;
			}
		};

		struct Decoder : public std::unary_function<const unsigned char, char> {
			const std::string decoding;
			Decoder(const std::string decoding) : decoding(decoding) {}
			char operator()(const unsigned char c) const {
				return decoding[c];
			}
		};
	
	};
	
};

#endif // __CHARMANIP_HH__
