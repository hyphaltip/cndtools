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
#include <iomanip>

#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;

#include "AlignmentReader.hh"
#include "LikelihoodCalculator.hh"

const std::string USAGE = "";

const std::string DESCRIPTION =
	"Given a multiple alignment and a selected sequence from the multiple "
	"alignment to be a putative recombinant, calculates the log likelihood of "
	"the data given the model at each parameter setting in POINTS";

std::istream& operator>>(std::istream& stream, std::vector<double>& p) {
	for (size_t i = 0; i < p.size(); ++i) {
		stream >> p[i];
	}
	return stream;
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	size_t recombinant_num = 0;
	std::string recombinant_name = "";
	std::string alignment_filename;
	std::string points_filename;
	
	// Set up option parser
	util::options::Parser parser(USAGE, DESCRIPTION);
	parser.addStoreOpt(0, "recombinant",
					   "Number (0-based) of sequence to treat as recombinant",
					   recombinant_num, "INTEGER");
	parser.addStoreOpt(0, "recombinant-name",
					   "Name of sequence to treat as recombinant",
					   recombinant_name, "STRING");
	parser.addStoreArg("alignmentFile",
					   "Alignment file",
					   alignment_filename);
	parser.addStoreArg("points",
					   "File containing points at which to calculate "
					   "the likelihood",
					   points_filename);
	parser.parse(argv, argv + argc);

	try {
		// Read in alignment
		InputFileStream alignment_file(alignment_filename);
		AlignmentReader reader(alignment_file,
							   recombinant_name, recombinant_num);

		InputFileStream points_file(points_filename);
		LikelihoodCalculator lc;

		std::vector<double> point(2);
		std::cout << std::setprecision(20);
		while (points_file >> point) {
			std::cout << lc.likelihood(reader.getAlignment(),
									   reader.getRecombinantNum(),
									   point[0],
									   point[1])
					  << '\n';
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
