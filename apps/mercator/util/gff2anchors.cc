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
#include <stdexcept>

#include "bio/gff/GFFRecord.hh"
#include "bio/gff/GFFInputStream.hh"
#include "util/options.hh"
#include "util/parser.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Values to be (optionally) specified on the command line
	int minLength = 5;
	
	// Parse options
	util::options::Parser parser("< gffInput > anchorOutput",
								 "Output genomic anchors based on CDS "
								 "GFF records");
	parser.addStoreOpt(0, "min-length",
					   "Minimum anchor length (# of amino acids)",
					   minLength, "LENGTH");
	parser.parse(argv, argv + argc);

	try {
		// Construct GFF stream for fast reading
		bio::gff::GFFInputStream strm(std::cin);

		// Write each GFF record to standard output
		size_t anchorNum = 0;
		bio::gff::GFFRecord rec;
		while (strm >> rec) {
			if (rec.getFeature() == "CDS" and rec.hasFrame()) {
				bio::genome::Position start, end;
				if (rec.getStrand() == '+') {
					start = rec.getStart() - 1 + rec.getFrame();
					end = start + ((rec.getEnd() - start) / 3) * 3;
				} else {
					end = rec.getEnd() - rec.getFrame();
					start = end - ((end - (rec.getStart() - 1)) / 3) * 3;
				}

				// Continue if length will be less than one codon's worth
				if (end < start or (end - start) < (minLength * 3)) {
					continue;
				}

				// Output anchor line
				std::cout << anchorNum << '\t'
						  << rec.getSeqname() << '\t'
						  << rec.getStrand() << '\t'
						  << start << '\t'
						  << end << '\t'
						  << 1 << '\n';

				++anchorNum;
			}
		}
	}
	// Report any formatting errors
	catch (const util::parser::FormatError& e) {
		std::cerr << "Error: " << e.getProblem() << '\n'
				  << e.getLine() << '\n';
		return EXIT_FAILURE;
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
