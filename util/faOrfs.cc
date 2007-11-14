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
#include <string>

#include "bio/formats/fasta/InputStream.hh"
#include "bio/gff/GFFRecord.hh"
#include "bio/alphabet/AmbiguousNucleotide.hh"
#include "bio/translation/Translator.hh"
#include "util/options.hh"
using bio::translation::Translator;

void findOrfs(std::string& seq,
			  const std::string& seqname,
			  const Translator& translator,
			  const unsigned int minLength,
			  const bool mask,
			  std::ostream& strm) {
	if (mask) {
		bio::alphabet::Nucleotide::hardMaskInPlace(seq);
	}

	bio::gff::GFFRecord rec;
	rec.setSeqname(seqname);
	rec.setSource("faOrfs");
	rec.setFeature("CDS");
	rec.setFrame(0);
	
	rec.setStrand('+');
	for (unsigned int phase = 0; phase < 3; ++phase) {
		std::string protein = translator.translate(seq, phase);

		std::string::size_type start = protein.find_first_not_of("X*");
		while (start != std::string::npos) {
			std::string::size_type end = protein.find_first_of("X*", start + 1);
			unsigned int length = 
				end == std::string::npos ?
				(protein.length() - start) * 3:
				(end - start) * 3;
			if (length >= minLength) {
				rec.setStart(start * 3 + phase + 1);
				rec.setEnd(rec.getStart() + (length - 1));
				strm << rec;
			}
			if (end == std::string::npos) {
				break;
			} else {
				start = protein.find_first_not_of("X*", end + 1);
			}
		}
	}

	bio::alphabet::AmbiguousDNA.reverseComplementInPlace(seq);
	
	rec.setStrand('-');
	for (unsigned int phase = 0; phase < 3; ++phase) {
		std::string protein = translator.translate(seq, phase);
		std::string::size_type start = protein.find_first_not_of("X*");	
		while (start != std::string::npos) {
			std::string::size_type end = protein.find_first_of("X*", start + 1);
			unsigned int length = 
				end == std::string::npos ?
				(protein.length() - start) * 3:
				(end - start) * 3;
			if (length >= minLength) {
				rec.setEnd(seq.length() - (start * 3 + phase));
				rec.setStart(rec.getEnd() - (length - 1));
				strm << rec;
			}
			if (end == std::string::npos) {
				break;
			} else {
				start = protein.find_first_not_of("X*", end + 1);
			}
		}
	}
	
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	size_t tableNum = 1;
	bool mask = false;
	size_t minLength = 150;

	util::options::Parser parser("< fastaInput > gffOutput",
								 "Output all potential ORFs in FASTA input");
	parser.addStoreOpt('t', "table",
					   "Number of the translation table to use for translation",
					   tableNum, "NUM");
	parser.addStoreTrueOpt('m', "mask",
						   "do not predict ORFs in masked regions",
						   mask);
	parser.addStoreOpt('l', "minlength", 
					   "only output ORFs of length at least LENGTH (bp)",
					   minLength, "LENGTH");
	parser.parse(argv, argv + argc);

	try {
		Translator translator(tableNum);

		// Construct FASTA stream for fast reading
		bio::formats::fasta::InputStream fastaStream(std::cin);
		
		// Step through records, output ORFs for each
		bio::formats::fasta::Record rec;
		while (fastaStream >> rec) {
			findOrfs(rec.sequence, rec.title,
					 translator, minLength, mask, std::cout);
		}
	
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		exit(EXIT_FAILURE);
	}
	
	return EXIT_SUCCESS;
}
