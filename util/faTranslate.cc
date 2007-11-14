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
#include "bio/translation/Translator.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "util/options.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	size_t frame = 0;
	bool toStop = false;
	size_t tableNum = 1;

	// Parse options
	util::options::Parser parser("< fastaInput",
								 "Translate DNA FASTA records to protein");
	parser.addStoreOpt('f', "frame",
					   "Frame (0, 1, or 2) in which to translate each FASTA "
					   "record",
					   frame, "NUM");
	parser.addStoreOpt('t', "table",
					   "Number of the translation table to use for translation",
					   tableNum, "NUM");
	parser.addStoreTrueOpt('s', "to-stop",
						   "Translate each DNA record to the first stop codon",
						   toStop);
	parser.parse(argv, argv + argc);

	try {
	
		// Validate options
		if (frame < 0 or frame > 2) {
			throw std::runtime_error("Frame must be 0, 1, or 2");
		}

		bio::translation::Translator translator(tableNum);

		// Create FASTA stream for fast reading
		bio::formats::fasta::InputStream fastaInStream(std::cin);
		bio::formats::fasta::OutputStream fastaOutStream(std::cout);

		// Step through records, translating and outputting
		bio::formats::fasta::Record rec;
		while (fastaInStream >> rec) {
			rec.sequence = translator.translate(rec.sequence,
												frame,
												toStop);
			fastaOutStream << rec;
		}
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
