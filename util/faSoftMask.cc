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

#include <iostream>  // for I/O
#include <algorithm> // for transform()
#include <cctype>    // for tolower()

#include "bio/formats/fasta.hh"
#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;
using namespace bio::formats;

struct SoftMasker {
	char operator()(const char c, const char masker) const {
		return std::toupper(masker) == 'N' ? std::tolower(c) : c;
	}
};

std::string softMask(const std::string& unmasked,
					 const std::string& hardmasked) {
	std::string softmasked;
	softmasked.reserve(unmasked.size());
	std::transform(unmasked.begin(), unmasked.end(),
				   hardmasked.begin(),
				   std::back_inserter(softmasked),
				   SoftMasker());
	return softmasked;
}

fasta::Record softMaskFASTA(fasta::Record& unmasked,
										  fasta::Record& hardmasked) {
	return fasta::Record(unmasked.title,
									   softMask(unmasked.sequence,
												hardmasked.sequence));
}

void softMaskFASTAStream(fasta::InputStream& unmasked,
						 fasta::InputStream& hardmasked,
						 fasta::OutputStream& softmasked) {
	while (true) {
		fasta::Record unmaskedRec, hardmaskedRec;

		unmasked >> unmaskedRec;
		hardmasked >> hardmaskedRec;

		// Check for reaching end of file in one or both files
		if (not unmasked or not hardmasked) {
			if (unmasked or hardmasked) {
				throw std::runtime_error("FASTA files have different "
										 "numbers of records");
			}
			return;
		}

		// Check that the two records have the same title
		if (unmaskedRec.title != hardmaskedRec.title) {
			throw std::runtime_error("Different records reached: " +
									 unmaskedRec.title +
									 " and " +
									 hardmaskedRec.title);
		}

		// Check that the two records have the same length
		if (unmaskedRec.sequence.length() != hardmaskedRec.sequence.length()) {
			throw std::runtime_error("Record titled " +
									 unmaskedRec.title +
									 "has different lengths in the two files");
		}
		
		// Soft mask using the hard masked sequence and write out the record
		softmasked << softMaskFASTA(unmaskedRec, hardmaskedRec);
	}
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Arguments
	std::string unmaskedFilename;
	std::string hardmaskedFilename;
	
	// Parse command line
	util::options::Parser parser("> softmaskedFASTAOutput",
								 "Produce a softmasked FASTA file from "
								 "hardmasked and unmasked files");
	parser.addStoreArg("unmaskedFile",
					   "FASTA file containing unmasked sequence",
					   unmaskedFilename);
	parser.addStoreArg("hardmaskedFile",
					   "FASTA file containing hardmasked sequence",
					   hardmaskedFilename);
	parser.parse(argv, argv + argc);

	try {
		InputFileStream unmaskedFile(unmaskedFilename);
		fasta::InputStream unmaskedStream(unmaskedFile);
		
		InputFileStream hardmaskedFile(hardmaskedFilename);
		fasta::InputStream hardmaskedStream(hardmaskedFile);	
		
		fasta::OutputStream outputStream(std::cout);
		
		softMaskFASTAStream(unmaskedStream, hardmaskedStream, outputStream);

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
