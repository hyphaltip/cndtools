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
#include <map>
#include <fstream>

#include "bio/formats/fasta/InputStream.hh"
#include "util/options.hh"
#include "util/io/line/InputStream.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	std::string coordsFilename;
	
	util::options::Parser parser("< fastaInput > alignOutput",
								 "Convert FASTA alignment to ALIGN format");
	parser.addStoreOpt(0, "coords",
					   "File containing sequence coordinates for the records "
					   "in the fast input",
					   coordsFilename, "FILE");
	parser.parse(argv, argv + argc);

	typedef std::map<std::string, std::string> CoordMap;
	CoordMap coordMap;
	if (!coordsFilename.empty()) {
		std::ifstream coordsFile(coordsFilename.c_str());
		util::io::line::InputStream lineStream(coordsFile);
		std::string line;
		while (lineStream >> line) {
			std::string::size_type tabPos = line.find_first_of('\t');
			if (tabPos == std::string::npos) {
				continue;
			}
			coordMap[line.substr(0, tabPos)] = line.substr(tabPos + 1);
		}
	}
	
	// Construct FASTA stream for fast reading
	bio::formats::fasta::InputStream fastaStream(std::cin);
	
	bio::formats::fasta::Record rec;
	while (fastaStream >> rec) {
		CoordMap::const_iterator it = coordMap.find(rec.title);
		std::cout << (it == coordMap.end() ? rec.title : it->second) << '\t';

		// Check for empty sequence
		if (rec.sequence.length() == 0) {
			std::cout << '\n';
			continue;
		}

		std::string::size_type seqStart = 0; // Start of ungapped sequence
		std::string::size_type seqEnd = 0; // End of ungapped sequence
		std::string::size_type gapEnd = 0; // End of gap after sequence
		while (true) {
			// Find next sequence and gap indices
			seqStart = gapEnd;
			seqEnd = rec.sequence.find_first_of("-", seqStart);
			gapEnd = rec.sequence.find_first_not_of("-", seqEnd);

			// Adjust coordinates at the end of the sequence
			if (seqEnd == std::string::npos) {
				seqEnd = rec.sequence.length();
			}
			if (gapEnd == std::string::npos) {
				gapEnd = rec.sequence.length();
			}

			// Output sequence length and gap length
			std::cout << seqEnd - seqStart << ',' << gapEnd - seqEnd;

			// Stop if the last gap ends at the end of the sequence
			if (gapEnd == rec.sequence.length()) {
				std::cout << '\n';
				break;
			} else {
				std::cout << ',';
			}
		}
	}

	return EXIT_SUCCESS;
}
