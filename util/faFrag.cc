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

#include <iostream>
#include <vector>

#include "bio/formats/fasta/InputStream.hh"
#include "bio/formats/fasta/OutputStream.hh"
#include "util/string.hh"
#include "util/options.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	std::vector<int> coords;
	util::options::Parser parser("start [end] < in.fa > out.fa", "");
	parser.addAppendArg("", "", coords, 1, 2);
	parser.parse(argv, argv + argc);

	try {
		// Construct FASTA stream for fast reading and read in first record
		bio::formats::fasta::InputStream fastaStream(std::cin);
		bio::formats::fasta::Record rec;
		fastaStream >> rec;
		
		if (!fastaStream) {
			throw std::runtime_error("Error reading fasta record from input");
		}

		int start, end;
		
		// Get start and end coordinates
		start = coords[0];
		end = (coords.size() == 2 ? coords[1] : rec.sequence.length());
		
		// Adjust negative coordinates relative to sequence end
		start = (start < 0 ? start + rec.sequence.length() : start);
		end = (end < 0 ? end + rec.sequence.length() : end);
		
		// Check that coordinates are valid
		if (start < 0 ||
			end < 0 ||
			static_cast<unsigned int>(start) > rec.sequence.length() ||
			static_cast<unsigned int>(end) > rec.sequence.length() ||
			start > end) {
			throw std::runtime_error("Invalid coordinates");
		}
		
		// Change the title
		rec.title = rec.title + ":" +
			util::string::toString(start) + "-" + util::string::toString(end);
		// Extract the substring
		rec.sequence = rec.sequence.substr(start, end - start);
		
		// Write the modified FASTA record to stdout
		bio::formats::fasta::OutputStream fastaOutStream(std::cout);
		fastaOutStream << rec;
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
