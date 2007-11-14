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

#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;

#include "AlignmentReader.hh"
#include "BestStateSetLister.hh"

const std::string USAGE = "";

const std::string DESCRIPTION =
"Given a multiple alignment, a selected sequence from the multiple "
"alignment to be a putative recombinant, and parameter settings, "
"returns the list of sequence subsets making up all optimal paths according "
"to the recombination model.";

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	size_t recombinant_num = 0;
	std::string recombinant_name = "";
	int match = 0;
	int mismatch = -2;
	int recombination = -8;
	int norecombination = 0;
	std::string alignment_filename;
	
	// Set up option parser
	util::options::Parser parser(USAGE, DESCRIPTION);
	parser.addStoreOpt(0, "recombinant",
					   "Number (0-based) of sequence to treat as recombinant",
					   recombinant_num, "INTEGER");
	parser.addStoreOpt(0, "recombinant-name",
					   "Name of sequence to treat as recombinant",
					   recombinant_name, "STRING");
	parser.addStoreOpt('m', "match", "Match score", match, "INTEGER");
	parser.addStoreOpt('x', "mismatch", "Mismatch score", mismatch, "INTEGER");
	parser.addStoreOpt('r', "recombination", "Recombination score",
					   recombination, "INTEGER");
	parser.addStoreOpt('n', "norecombination", "Norecombination score",
					   norecombination, "INTEGER");
	parser.addStoreArg("alignmentFile",
					   "Alignment file",
					   alignment_filename);
	parser.parse(argv, argv + argc);

	try {
		InputFileStream alignment_file(alignment_filename);
		AlignmentReader reader(alignment_file,
							   recombinant_name, recombinant_num);

		BestStateSetLister lister(match, mismatch,
								  recombination, norecombination);
		StateSetList bestStateSetList;
		lister.getBestStateSetList(reader.getAlignment(),
								   reader.getRecombinantNum(),
								   bestStateSetList);

		for (size_t i = 0; i < bestStateSetList.size(); ++i) {
			StateSet ss = bestStateSetList[i];
			size_t first = ss.find_first();
			for (size_t state = first; state < ss.size(); ++state) {
				if (not ss[state]) { continue; }
				if (state != first) { std::cout << ' '; }
				std::cout << state;
			}
			std::cout << '\n';
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
