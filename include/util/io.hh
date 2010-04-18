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

#ifndef __UTIL_IO_HH__
#define __UTIL_IO_HH__

#include <istream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdint.h>

#define IO_BUFFER_SIZE (4 * 1024)

namespace util { namespace io {
		
		void copy_stream(std::istream& in, std::ostream& out);
		
		std::string readStream(std::istream& stream);
		std::string readStream(FILE* stream);
		
		namespace binary {
			
			inline void swap_bytes(int8_t& x) { return; }

			inline void swap_bytes(uint8_t& x) { return; }

			inline void swap_bytes(int16_t& x) {		
				x = (((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8));
			}

			inline void swap_bytes(uint16_t& x) {		
				x = (((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8));
			}

			inline void swap_bytes(int32_t& x) {
				x = (((x & 0xFF000000L) >> 24) |               
					 ((x & 0x00FF0000L) >>  8) |              
					 ((x & 0x0000FF00L) <<  8) |              
					 ((x & 0x000000FFL) << 24));
			}

			inline void swap_bytes(uint32_t& x) {
				x = (((x & 0xFF000000L) >> 24) |               
					 ((x & 0x00FF0000L) >>  8) |              
					 ((x & 0x0000FF00L) <<  8) |              
					 ((x & 0x000000FFL) << 24));
			}

			inline void swap_bytes(int64_t& x) {
				union {
					int64_t ll;
					int32_t l[2];
				} w, r;
				w.ll = x;
				r.l[0] = w.l[1];
				r.l[1] = w.l[0];
				swap_bytes(r.l[0]);
				swap_bytes(r.l[1]);
				x = r.ll;
			}

			inline void swap_bytes(uint64_t& x) {
				union {
					uint64_t ll;
					uint32_t l[2];
				} w, r;
				w.ll = x;
				r.l[0] = w.l[1];
				r.l[1] = w.l[0];
				swap_bytes(r.l[0]);
				swap_bytes(r.l[1]);
				x = r.ll;
			}

			inline bool write(std::ostream& stream, const char* val, uint32_t len) {
				return stream.write(val, len);
			}
			
			inline bool write(FILE* stream, const char* val, uint32_t len) {
				return (fwrite(val, sizeof(char), len, stream) == len);
			}

			template<typename T, typename S>
			inline bool write(S& stream, T val) {
#if BYTE_ORDER == BIG_ENDIAN
		  	    swap_bytes(val);
#elif BYTE_ORDER != LITTLE_ENDIAN
#error "Byte order macros not defined"
#endif
				return write(stream, reinterpret_cast<const char*>(&val), sizeof(T));
			}

			inline bool read(std::istream& stream, char* val, uint32_t len) {
				return stream.read(val, len);
			}

			inline bool read(FILE* stream, char* val, uint32_t len) {
				return (fread(val, sizeof(char), len, stream) == len);
			}

			template<typename S, typename T>
			inline bool read(S& stream, T& val) {
				bool result = read(stream, reinterpret_cast<char*>(&val), sizeof(T));
#if BYTE_ORDER == BIG_ENDIAN
		  	    swap_bytes(val);
#elif BYTE_ORDER != LITTLE_ENDIAN
#error "Byte order macros not defined"
#endif
				return result;
			}
		
			template<typename S>
			inline bool write(S& stream, const std::string& val) {
				uint32_t len = val.length();
				return write(stream, len) && write(stream, val.data(), len);
			}
	
			template<typename S>
			inline bool read(S& stream, std::string& val) {
				uint32_t len = 0;
				if (read(stream, len)) {
					std::vector<char> buffer(len);
					if (read(stream, &buffer[0], len)) {
						val.assign(&buffer[0], len);
						return true;
					} else {
						return false;
					}
				} else {
					return false;
				}
			}
		
			template<typename T>
			inline uint32_t size(const T&) { return sizeof(T); }

			template<>
			inline uint32_t size<std::string>(const std::string& val) {
				return sizeof(uint32_t) + val.length();
			}

		};

	}  }

#endif // __UTIL_IO_HH__
