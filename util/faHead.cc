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

#include "bio/formats/fasta.hh"
#include "util/string.hh"
#include "util/options.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	int numRecs = 10;

	util::options::Parser parser("< fastaInput",
								 "Output the first N records of fasta input");
	parser.addStoreOpt('n', "numrecs",
					   "number of records to output",
					   numRecs, "NUM");
	parser.parse(argv, argv + argc);

	// Construct FASTA stream for fast reading
	bio::formats::fasta::InputStream fastaInStream(std::cin);
	bio::formats::fasta::OutputStream fastaOutStream(std::cout);

	// Step through records
	bio::formats::fasta::Record rec;
	int counter = 0;
    while (counter < numRecs && fastaInStream >> rec) {
		fastaOutStream << rec;
		++counter;
	}

	return EXIT_SUCCESS;
}
