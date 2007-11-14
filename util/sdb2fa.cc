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
#include <string>

#include "bio/formats/fasta.hh"
#include "bio/sdb.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "bio/alphabet/AminoAcid.hh"
#include "util/string.hh"
#include "util/options.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool number_title = false;
	bool unmask = false;
	bool hardmask = false;
	bool protein = false;
	std::string dbFilename;

	util::options::Parser parser("", "");
	parser.addStoreTrueOpt('t', "number-titles", 
						   "use record numbers as FASTA title in output",
						   number_title);
	parser.addStoreTrueOpt('u', "unmask",
						   "output all sequences in upper case (lose softmask formatting)",
						   unmask);
	parser.addStoreTrueOpt('h', "hardmask",
						   "hardmask the sequences.  Masked residues will be replaced by 'N' for DNA (default) and 'X' for protein",
						   hardmask);
	parser.addStoreTrueOpt('p', "protein",
						   "treat sequences as protein.  Only used in combination with the hardmask option",
						   protein);
	parser.addStoreArg("dbFile", "", dbFilename);
	parser.parse(argv, argv + argc);
	
	if (unmask && hardmask) {
		std::cerr << "Error: "
				  << "Cannot specify both unmask and hardmask options\n";
		exit(EXIT_FAILURE);
	}

	bio::SDB::DB db;
	
	try {
		db.open(dbFilename, true);
	} catch (std::runtime_error& err) {
		std::cerr << "Error: Could not open database: " << err.what() << '\n';
		exit(EXIT_FAILURE);
	}

	bio::formats::fasta::OutputStream fastaOutStream(std::cout);
	bio::formats::fasta::Record rec;
	for (bio::SDB::DB::Iterator it = db.begin(); it != db.end(); ++it) {
		rec.title = number_title ?
			util::string::toString(it->getRecNum()) : it->getTitle();
		rec.sequence = it->getSeq();
		if (unmask) {
			if (protein) {
				bio::alphabet::AminoAcid::unMaskInPlace(rec.sequence);
			} else {
				bio::alphabet::Nucleotide::unMaskInPlace(rec.sequence);
			}
		} else if (hardmask) {
			if (protein) {
				bio::alphabet::AminoAcid::hardMaskInPlace(rec.sequence);
			} else {
					bio::alphabet::Nucleotide::hardMaskInPlace(rec.sequence);
			}
		}
		fastaOutStream << rec;
	}

	return EXIT_SUCCESS;
}
