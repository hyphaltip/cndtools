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
#include <iterator>
#include <vector>
#include <algorithm>

#include "bio/formats/fasta.hh"
#include "util/stl.hh"
#include "util/options.hh"

// Functor for comparing FASTA records based on their titles
struct FASTATitleComparer {
	bool operator()(const bio::formats::fasta::Record& rec1,
					const bio::formats::fasta::Record& rec2) const {
		return rec1.title < rec2.title;
	}
};

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Parse command line
	util::options::Parser parser("< fastaInput", "Sort FASTA records by title");
	parser.parse(argv, argv + argc);

	try {
		// Construct FASTA stream for fast reading
		bio::formats::fasta::InputStream fastaInStream(std::cin);
		
		// Read FASTA records into a vector
		std::vector<bio::formats::fasta::Record> recs;
		bio::formats::fasta::Record rec;
		while (fastaInStream >> rec) {
			recs.push_back(rec);
		}
		
		// Sort the records based on their titles
		std::sort(recs.begin(), recs.end(), FASTATitleComparer());
		
		// Write the sorted records to stdout
		bio::formats::fasta::OutputStream fastaOutStream(std::cout);
		typedef std::vector<bio::formats::fasta::Record>::const_iterator Iterator;
		for (Iterator it = recs.begin(); it != recs.end(); ++it) {
			fastaOutStream << *it;
		}
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
		
	return EXIT_SUCCESS;
}
