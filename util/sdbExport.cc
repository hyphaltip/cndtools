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
#include <stdexcept>

#include "bio/formats/fasta.hh"
#include "bio/sdb.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "bio/alphabet/AminoAcid.hh"
#include "util/string.hh"
#include "util/options.hh"
#include "util/io/line/InputStream.hh"
using namespace bio;
using namespace bio::formats;
using namespace util::io;
using util::string::toString;

typedef std::vector<std::string> Coords;

util::string::Converter<int> toInteger;

// Get the record from the database according to name or number
void getRecord(SDB::Record& rec,
			   SDB::DB& db,
			   Coords& coords,
			   bool number) {
	if (number) {
		size_t recNum = toInteger(coords[0]);
		db.getRec(recNum, rec);
	} else {
		db.getRec(coords[0], rec);
	}
}

std::string toString(Coords& coords) {
	return util::string::join(coords.begin(), coords.end(), " ");
}

void makeFASTARecord(fasta::Record& fastaRec,
					 SDB::Record& sdbRec,
					 Coords& coords,
					 const std::string& name,
					 bool number_title,
					 bool range) {
	// Initialize the coordinates
	size_t start = 0;
	size_t end = sdbRec.getLength();
	
	// Adjust coordinates
	if (coords.size() > 1) {
		int userStart = toInteger(coords[1]);
		start = userStart >= 0 ? userStart : userStart + sdbRec.getLength();
	}
	if (coords.size() > 2) {
		int userEnd = toInteger(coords[2]);
		end = userEnd >= 0 ? userEnd : userEnd + sdbRec.getLength();
	}
	
	// Check the validity of the coordinates
	if (start > end or start > sdbRec.getLength() or end > sdbRec.getLength()) {
		throw std::runtime_error("Invalid coordinates: " +
								 toString(coords));
	}

	// Retrieve sequence for FASTA Record
	fastaRec.sequence = sdbRec.getSeq(start, end);
	
	// Determine the title of the FASTA record
	if (not name.empty()) {
		fastaRec.title = name;
	} else if (number_title) {
		fastaRec.title = toString(sdbRec.getRecNum());
	} else {
		fastaRec.title = sdbRec.getTitle();
	}
	if (range) {
		fastaRec.title += ":" + toString(start) + "-" + toString(end);
		if (coords.size() > 3) {
			fastaRec.title += coords[3];
		}
	}
}

void processFASTARecord(fasta::Record& fastaRec,
						Coords& coords,
						bool unmask,
						bool hardmask,
						bool protein) {
	// Reverse complement sequence if strand is "-"
	if (coords.size() == 4) {
		if (protein) {
			throw std::runtime_error("Specified strand for protein");
		}
		if (coords[3] == "-") {
			alphabet::AmbiguousDNA.reverseComplementInPlace(fastaRec.sequence);
		} else if (coords[3] != "+") {
			throw std::runtime_error("Invalid strand: must be + or -");
		}
	}

	// Mask sequence as specified
	if (unmask) {
		if (protein) {
			alphabet::AminoAcid::unMaskInPlace(fastaRec.sequence);
		} else {
			alphabet::Nucleotide::unMaskInPlace(fastaRec.sequence);
		}
	} else if (hardmask) {
		if (protein) {
			alphabet::AminoAcid::hardMaskInPlace(fastaRec.sequence);
		} else {
			alphabet::Nucleotide::hardMaskInPlace(fastaRec.sequence);
		}
	}
}

line::InputStream& operator>>(line::InputStream& stream, Coords& coords) {
	static std::string line;
	if (stream >> line) {
		coords.clear();
		util::string::split(line, std::back_inserter(coords));
		if (coords.size() == 0 or coords.size() > 4) {
			throw std::runtime_error("Invalid coordinates: " +
									 toString(coords));
		}
	}
	return stream;
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::string name;
	bool number = false;
	bool number_title = false;
	bool range = false;
	bool unmask = false;
	bool hardmask = false;
	bool protein = false;
	std::string dbFilename;
	std::vector<std::string> coords;
	
	util::options::Parser parser("[recName [start [end [strand]]]] "
								 "[< intervals] > FASTARecords", "");
	parser.addStoreOpt(0, "name",
					   "Use NAME for the title of the record",
					   name, "NAME");
	parser.addStoreTrueOpt('n', "number",
						   "recName specifies a record number instead of a "
						   "record title",
						   number);
	parser.addStoreTrueOpt('t', "number-titles",
						   "use record numbers as FASTA title in output",
						   number_title);
	parser.addStoreTrueOpt('r', "range",
						   "append the range 'start-end' to the FASTA title",
						   range);
	parser.addStoreTrueOpt('u', "unmask",
						   "output entire sequence in upper case (lose "
						   "softmask formatting)",
						   unmask);
	parser.addStoreTrueOpt('h', "hardmask",
						   "hardmask the sequence.  Masked residues will be "
						   "replaced by 'N' for DNA (default) and 'X' for "
						   "protein",
						   hardmask);
	parser.addStoreTrueOpt('p', "protein",
						   "treat sequence as protein.  Only used in "
						   "combination with the hardmask option",
						   protein);
	parser.addStoreArg("dbFile", "", dbFilename);
	parser.addAppendArg("", "", coords, 0, 4);
	parser.parse(argv, argv + argc);
	
	try {
		if (unmask and hardmask) {
			throw std::runtime_error("Cannot specify both unmask and "
									 "hardmask options");
		}
		
		SDB::DB db;
		db.open(dbFilename);

		fasta::OutputStream fastaOutStream(std::cout);
		
		SDB::Record rec;
		fasta::Record faRec;
		
		if (coords.empty()) {
			db.readIndex();
			util::io::line::InputStream lineStream(std::cin);
			while (lineStream >> coords) {
				getRecord(rec, db, coords, number);
				makeFASTARecord(faRec, rec, coords, name, number_title, range);
				processFASTARecord(faRec, coords, unmask, hardmask, protein);
				fastaOutStream << faRec;
			}
		} else {
			getRecord(rec, db, coords, number);
			makeFASTARecord(faRec, rec, coords, name, number_title, range);
			processFASTARecord(faRec, coords, unmask, hardmask, protein);
			fastaOutStream << faRec;
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
