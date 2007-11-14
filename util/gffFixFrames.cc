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
#include <fstream>

#include "bio/gff/GFFRecord.hh"
#include "bio/gff/GFFInputStream.hh"
#include "bio/sdb.hh"
#include "bio/formats/fasta/InputStream.hh"
#include "bio/translation/Table.hh"
#include "util/parser.hh"
#include "util/stl.hh"
#include "util/options.hh"

int main(int argc, const char** argv) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);
	
	// Values to be (optionally) specified on the command line
	std::string dbFilename;
	std::string protFilename;
	std::string badRecFilename;
	std::string attribute = "name";
	size_t tableNum = 1;

	util::options::Parser parser("< gffInput > gffOutput",
								 "Correct the frames of CDS records in GFF input");
	parser.addStoreOpt('t', "table",
					   "Number of the translation table to use for translation",
					   tableNum, "NUM");
	parser.addStoreOpt('i', "id",
					   "gff attribute to link CDS records with protein records",
					   attribute, "ATTRIBUTE");
	parser.addStoreOpt(0, "badrecs",
					   "file to output bad records",
					   badRecFilename, "FILENAME");
	parser.addStoreArg("proteins.fa",
					   "FASTA file containing the protein sequences of the "
					   "genes specified in the GFF input",
					   protFilename);
	parser.parse(argv, argv + argc);

	try {
		const bio::translation::Table* table =
			bio::translation::Table::getTable(tableNum);
		if (table == NULL) {
			throw std::runtime_error("Invalid translation table: " +
									 util::string::toString(tableNum));
		}
		
		// Open the database
		bio::SDB::DB db;
		db.open(dbFilename, true);
		
		// Attempt to open the protein file
		std::ifstream protFile(protFilename.c_str());
		if (!protFile) {
			throw std::runtime_error("Could not open protein file " +
									 protFilename);
		}
		
		// Attempt to open the bad record file
		std::ofstream* badFile = (badRecFilename.empty() ?
								  NULL :
								  new std::ofstream(badRecFilename.c_str()));
		if (badFile != NULL && !(*badFile)) {
			throw std::runtime_error("Could not open bad record file " +
									 badRecFilename);
		}
		
		// Construct a FASTA stream for fast reading
		bio::formats::fasta::InputStream fastaStream(protFile);
		
		// Read the protein fasta file into a hash table
		typedef util::stl::hash_map<std::string, bio::formats::fasta::Record> ProteinMap;
		ProteinMap proteins;
		bio::formats::fasta::Record faRec;
		while (fastaStream >> faRec) {
			proteins[faRec.title] = faRec;
		}
	
		// Construct GFF stream for fast reading
		bio::gff::GFFInputStream gffStream(std::cin);

		bio::gff::GFFRecord rec;
		bio::SDB::Record chromRec;
		while (gffStream >> rec) {
			bool bad = false;
			
			// Find sequence in database
			try {
				db.getRec(rec.getSeqname(), chromRec);
			} catch (bio::SDB::NotFoundError& err) {
				std::cerr << err.what() << '\n';
				bad = true;
			}

			if (bad) {
				// pass
			}
			// Check to make sure the coordinates are ok
			else if (rec.getStart() < 1 ||
					 rec.getEnd() > chromRec.getLength() ||
					 rec.getStart() > rec.getEnd()) {
				std::cerr << "Error: Invalid coordinates: " << rec;
				bad = true;
			}
			// Check start codon
			else if (rec.getFeature() == "start_codon") {
				std::string dna = chromRec.getSeq(rec.getStart() - 1,
												  rec.getEnd(),
												  rec.getStrand());
				if (!table->isStartCodon(dna)) {
					std::cerr << "Error: Bad start codon ("
							  << dna << "): " << rec;
					bad = true;
				}
			}
			// Check stop codon
			else if (rec.getFeature() == "stop_codon") {
				std::string dna = chromRec.getSeq(rec.getStart() - 1,
												  rec.getEnd(),
												  rec.getStrand());
				if (!table->isStopCodon(dna)) {
					std::cerr << "Error: Bad stop codon ("
							  << dna << "): " << rec;
					bad = true;
				}
			}
			// Check CDS
			else if (rec.getFeature() == "CDS") {
				// Find the protein sequence corresponding to this CDS
				ProteinMap::iterator pos = proteins.find(rec.getAttribute(attribute).values.front());
				if (pos == proteins.end()) {
					std::cerr << "Error: No protein with title "
							  << rec.getFeature() << '\n';
					bad = true;
				} else {
					std::string protein = pos->second.sequence;
					
					// Extract the sequence for this CDS
					std::string dna = chromRec.getSeq(rec.getStart() - 1,
													  rec.getEnd(),
													  rec.getStrand());
					
					// Translate according to the given frame
					std::string aas = table->translate(dna, rec.getFrame());
					// Check if this was the correct frame
					bad = (protein.find(aas) == std::string::npos);
					if (bad) {
						std::cerr << "Error: Incorrect frame for: " << rec;
						for (unsigned int f = 0; f < 3; ++f) {
							if (f == rec.getFrame()) {
								continue;
							}
							aas = table->translate(dna, f);
							if (protein.find(aas) != std::string::npos) {
								rec.setFrame(f);
								bad = false;
								break;
							}
						}
					}
				}
			}

			// Output record to bad record file
			if (bad && badFile) {
				(*badFile) << rec;
			}

			// Output (corrected) record
			std::cout << rec;
		}
	} catch (util::parser::FormatError& e) {
		std::cerr << "Error: " << e.getProblem() << '\n' << e.getLine() << '\n';
		return EXIT_FAILURE;
	} catch (std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
