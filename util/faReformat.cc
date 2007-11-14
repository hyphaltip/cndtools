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
using namespace bio::formats::fasta;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool stripTitle = false;
	bool removeGaps = false;
	bool removeEmpty = false;
	size_t width = Constants::DEFAULT_LINE_WIDTH;

	// Parse command line
	util::options::Parser parser("< fastaInput", "Reformat fasta input");
	parser.addStoreTrueOpt('s', "strip-title",
						   "truncate the title for each record at the "
						   "first whitespace character",
						   stripTitle);
	parser.addStoreTrueOpt('g', "remove-gaps",
						   "remove gaps from sequences",
						   removeGaps);
	parser.addStoreTrueOpt('e', "remove-empty",
						   "remove empty records",
						   removeEmpty);
	parser.addStoreOpt('w', "width",
					   "number of characters at which sequence is wrapped",
					   width, "NUM");
	parser.parse(argv, argv + argc);
	
	InputStream fastaInStream(std::cin);

	OutputStream fastaOutStream(std::cout);
	fastaOutStream.setLineWidth(width);
	
	// Step through records, possibly cleaning, and then output
	Record rec;
    while (fastaInStream >> rec) {
		if (stripTitle) {
			rec.title = util::string::firstWord(rec.title);
		}
		if (removeGaps) {
			rec.sequence.erase(std::remove(rec.sequence.begin(),
										   rec.sequence.end(),
										   '-'),
							   rec.sequence.end());
		}
		if (not removeEmpty or not rec.sequence.empty()) {
			fastaOutStream <<  rec;
		}
	}

	return EXIT_SUCCESS;
}
