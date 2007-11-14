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

#include <sstream>
#include <stdexcept>

#include "bio/formats/ucsc/net/InputStream.hh"

namespace bio { namespace formats { namespace ucsc { namespace net {

	class BadNetLineError : public std::runtime_error {
	public:
		BadNetLineError(const std::string& what,
						const std::string& line) :
			std::runtime_error("Bad net line (" + what + "):\n" + line) {
		}
	};
	
				
	InputStream& InputStream::operator>>(Record& r) {
		while (stream) {
			stream >> line;
			std::istringstream line_stream(line);

			// Read in the class of the line
			std::string rec_class;
			line_stream >> rec_class;
			
			if (not line_stream) { continue; }

			if (rec_class == "net") {
				line_stream >> tChrom >> tChromSize;
				if (not line_stream) {
					throw BadNetLineError("Invalid net record", line);
				}
				continue;
			} else if (rec_class == "fill") {
				r.recClass = Record::FILL;
			} else if (rec_class == "gap") {
				r.recClass = Record::GAP;
			} else {
				throw BadNetLineError("Unknown record class", line);
			}

			// Determine level from number of spaces
			r.level = line.find_first_not_of(' ');

			// Set fields defined by net record
			r.tChrom = tChrom;
			r.tChromSize = tChromSize;

			// Read in fixed fields
			line_stream >> r.tStart
						>> r.tSize
						>> r.qChrom
						>> r.orientation
						>> r.qStart
						>> r.qSize;
			if (not line_stream) {
				throw BadNetLineError("Invalid fill record", line);
			}

			// Read in tokens of attribute names and values
			std::string name, value;
			r.attributes.clear();
			while (line_stream >> name) {
				line_stream >> value;
				if (not line_stream) {
					throw BadNetLineError("No value for last attribute", line);
				}
				r.attributes[name] = value;
			}

			break;
		}

		return *this;
	}

} } } }
