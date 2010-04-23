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
#include <map>
#include <fstream>

#include "boost/tuple/tuple.hpp"
#include "boost/tokenizer.hpp"
using boost::tie;

#include "bio/formats/fasta/InputStream.hh"
#include "bio/formats/fasta/OutputStream.hh"
#include "bio/sdb.hh"
#include "bio/alphabet/AmbiguousNucleotide.hh"
#include "boost/unordered_map.hpp"
#include "util/io/line/InputStream.hh"
#include "util/io.hh"
#include "util/string.hh"
#include "util/options.hh"
#include "filesystem/Path.hh"

enum ExtractMethod {
	FASTA_FILE,
	SDB_FILE,
	SDB_DIR,
	PROG
};

class SeqExtractor {
public:
	SeqExtractor() {}
	virtual ~SeqExtractor() {}
	virtual void init(const std::string& arg) = 0;

	virtual std::string getSeq(const std::string& title) {
		throw std::runtime_error("Method does not support identifier format");
	}

	virtual std::string getSeq(const std::string& title,
							   const unsigned int start,
							   const unsigned int end,
							   const char strand) {
		throw std::runtime_error("Method does not support identifier format");
	}

	virtual std::string getSeq(const std::string& genome,
							   const std::string& chrom,
							   const unsigned int start,
							   const unsigned int end,
							   const char strand) {
		throw std::runtime_error("Method does not support identifier format");
	}
		
};

class SDBDirExtractor : public SeqExtractor {
private:
	filesystem::Path sdbDir;
public:
	void init(const std::string& arg) {
		sdbDir = arg;
	}

	std::string getSeq(const std::string& genome,
					   const std::string& chrom,
					   const unsigned int start,
					   const unsigned int end,
					   const char strand) {
		// Construct path to SDB file
		filesystem::Path sdbFilepath = sdbDir / (genome + ".sdb");
		
		// Open SDB file and extract seq
		bio::SDB::DB db;
		db.open(sdbFilepath);
		return db.getSeq(chrom, start, end, strand);
	}
};

class SDBExtractor : public SeqExtractor {
private:
	bio::SDB::DB db;
public:
	void init(const std::string& arg) {
		db.open(arg);
	}

	std::string getSeq(const std::string& title) {
		return db.getSeq(title);
	}

	std::string getSeq(const std::string& title,
					   const unsigned int start,
					   const unsigned int end,
					   const char strand) {
		return db.getSeq(title, start, end, strand);
	}
};

class ProgExtractor : public SeqExtractor {
private:
	std::string cmd;

	std::string runCommand(const std::string& fullcmd) {
		FILE* output = popen(fullcmd.c_str(), "r");
		std::istringstream strm(util::io::readStream(output));
		bio::formats::fasta::InputStream fastaStream(strm);
		bio::formats::fasta::Record rec;
		fastaStream >> rec;
		return rec.sequence;
	}

public:
	void init(const std::string& arg) {
		cmd = arg;
	}

	std::string getSeq(const std::string& title) {
		std::string fullcmd = cmd + " " + title;
		return runCommand(fullcmd);
	}
	
	std::string getSeq(const std::string& title,
					   const unsigned int start,
					   const unsigned int end,
					   const char strand) {
		std::string fullcmd = cmd + " " + title +
			" " + util::string::toString(start) +
			" " + util::string::toString(end) +
			" " + util::string::toString(strand);
		return runCommand(fullcmd);
	}
};

class FASTAExtractor : public SeqExtractor {
private:
    boost::unordered_map<std::string, std::string> seqMap;
public:
	void init(const std::string& arg) {
		std::ifstream fastaFile(arg.c_str());
		bio::formats::fasta::InputStream fastaStream(fastaFile);
		bio::formats::fasta::Record rec;
		while (fastaStream >> rec) {
			seqMap[rec.title] = rec.sequence;
		}
	}

	std::string getSeq(const std::string& title) {
		return seqMap[title];
	}

	std::string getSeq(const std::string& title,
					   const unsigned int start,
					   const unsigned int end,
					   const char strand) {
		std::string seq = seqMap[title].substr(start, end);
		if (strand == '-') {
			bio::alphabet::AmbiguousDNA.reverseComplementInPlace(seq);
		}
		return seq;
	}
};

std::pair<size_t, size_t> splitCoords(const std::string& coords) {
	util::string::Converter<size_t> toSizeT;
	std::string::size_type dashPos = coords.find_first_of('-');
	return std::make_pair(toSizeT(coords.substr(0, dashPos)),
						  toSizeT(coords.substr(dashPos + 1)));
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	ExtractMethod method = FASTA_FILE;
	std::string coordsFilename;
	std::string arg;

	util::options::Parser parser("fastaFile|sdbFile|sdbDir|command "
								 "< alignInput > fastaOutput", "");
	parser.addStoreConstOpt(0, "fasta",
							"extract sequence from multi-fasta file given as argument (DEFAULT). "
							"Sequence identifiers in align file must be of the form "
							"'recname:start-end:strand'",
							method, FASTA_FILE);
	parser.addStoreConstOpt(0, "sdb",
							"extract sequence from sdb file given as argument. "
							"Sequence identifiers in align file must be of the form "
							"'recname:start-end:strand'",
							method, SDB_FILE);
	parser.addStoreConstOpt(0, "sdbdir",
							"extract sequence from sdb files in directory given as argument. "
							"Sequence identifiers in align file must be of the form "
							"'sdbfileprefix:recname:start-end:strand'",
							method, SDB_DIR);
	parser.addStoreConstOpt(0, "prog",
							"extract sequence by running command given as argument",
							method, PROG);
	parser.addStoreOpt(0, "coords",
					   "File containing sequence coordinates for the records "
					   "in the fast input",
					   coordsFilename, "FILE");

	parser.addStoreArg("", "", arg);
	parser.parse(argv, argv + argc);

	try {
	
		typedef std::map<std::string, std::string> CoordMap;
		CoordMap coordMap;
		if (!coordsFilename.empty()) {
			std::ifstream coordsFile(coordsFilename.c_str());
			util::io::line::InputStream lineStream(coordsFile);
			std::string line;
			while (lineStream >> line) {
				std::string::size_type tabPos = line.find_first_of('\t');
				if (tabPos == std::string::npos) {
					continue;
				}
				coordMap[line.substr(tabPos + 1)] = line.substr(0, tabPos);
			}
		}
	
		SeqExtractor* extractor = NULL;
		switch (method) {
		case FASTA_FILE:
			extractor = new FASTAExtractor();
			break;
		case SDB_FILE:
			extractor = new SDBExtractor();
			break;
		case SDB_DIR:
			extractor = new SDBDirExtractor();
			break;
		case PROG:
			extractor = new ProgExtractor();
			break;
		}

		extractor->init(arg);
	
		util::io::line::InputStream lineStream(std::cin);
		bio::formats::fasta::OutputStream fastaOutStream(std::cout);
		std::string line;
		while (lineStream >> line) {
			std::string::size_type tabPos = line.find_first_of('\t');

			std::string nameString = line.substr(0, tabPos);
			std::string intervalString = line.substr(tabPos + 1);

			std::vector<size_t> intervals;
			typedef boost::tokenizer<boost::char_separator<char> > CharTok;
			CharTok splitter(intervalString, boost::char_separator<char>(","));
			std::transform(splitter.begin(), splitter.end(),
						   std::back_inserter(intervals),
						   util::string::Converter<size_t>());
		
			std::vector<std::string> nameTokens;
			util::string::split(nameString, std::back_inserter(nameTokens), ":");

			std::string seq;
			if (nameTokens.size() == 3) {
				size_t start, end;
				tie(start, end) = splitCoords(nameTokens[1]);
				seq = extractor->getSeq(nameTokens[0], start, end, nameTokens[2][0]);
			} else if (nameTokens.size() == 4) {
				size_t start, end;
				tie(start, end) = splitCoords(nameTokens[2]);
				seq = extractor->getSeq(nameTokens[0], nameTokens[1], start, end, nameTokens[3][0]);
			} else {
				seq = extractor->getSeq(nameString);
			}

			// Count gap lengths
			size_t totalGaps = 0;
			for (size_t i = 1; i < intervals.size(); i += 2) {
				totalGaps += intervals[i];
			}

			bio::formats::fasta::Record rec;

			CoordMap::const_iterator it = coordMap.find(nameString);
			rec.title = (it == coordMap.end() ? nameString : it->second);
		
			rec.sequence.reserve(seq.size() + totalGaps);

			std::string::const_iterator seqStart = seq.begin();

			for (size_t i = 0; i < intervals.size(); i += 2) {
				std::string::const_iterator seqEnd = seqStart + intervals[i];
				rec.sequence.append(seqStart, seqEnd);
				rec.sequence.append(intervals[i + 1], '-');
				seqStart = seqEnd;
			}
			
			fastaOutStream << rec;
		}
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
