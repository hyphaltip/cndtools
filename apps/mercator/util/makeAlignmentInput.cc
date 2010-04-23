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
#include <vector>
#include <map>
#include <algorithm>

#include "bio/sdb.hh"
#include "bio/genome/BasicInterval.hh"
#include "bio/formats/fasta.hh"
#include "bio/homologymap/Map.hh"
#include "bio/phylogenetic/Tree.hh"
#include "bio/formats/newick.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "util/io.hh"
#include "util/io/line/InputStream.hh"
#include "filesystem.hh"
using namespace bio;
using namespace filesystem;
using namespace bio::formats;
using util::string::toString;
using boost::unordered_map;
using boost::unordered_set;
using namespace util::io;

struct Anchor {
	std::string id;
	genome::BasicInterval interval;
	size_t isCoding;
};

typedef std::vector<Anchor*> AnchorList;
typedef std::map<std::string, Anchor*> AnchorMap;
typedef AnchorList Clique;
typedef std::vector<Clique> Run;
typedef unordered_map<size_t, Run> RunMap;
typedef std::vector<AnchorMap> GenomeAnchorMaps;

struct SimpleConstraint {
	std::string genome1;
	size_t start1;
	size_t end1;
	std::string genome2;
	size_t start2;
	size_t end2;
};

struct Constraint {
	size_t segNum;
	size_t isCoding;
	std::string genome1;
	genome::BasicInterval interval1;
	std::string genome2;
	genome::BasicInterval interval2;
};

std::ostream& operator<<(std::ostream& stream, const SimpleConstraint& c) {
	return stream << c.genome1 << ' '
				  << c.start1 << ' '
				  << c.end1 << ' '
				  << c.genome2 << ' '
				  << c.start2 << ' '
				  << c.end2 << '\n';
}

std::istream& operator>>(std::istream& stream, genome::BasicInterval& i) {
	std::string chrom;
	genome::Position start, end;
	genome::Strand strand;
	stream >> chrom >> start >> end >> strand;
	i.setChrom(chrom);
	i.setStart(start);
	i.setEnd(end);
	i.setStrand(strand);
	return stream;
}

std::istream& operator>>(std::istream& stream, Anchor& a) {
	return stream >> a.id >> a.interval >> a.isCoding;
}

std::istream& operator>>(std::istream& stream, Constraint& c) {
	return stream >> c.segNum
				  >> c.isCoding
				  >> c.genome1 >> c.interval1
				  >> c.genome2 >> c.interval2;
}

std::ostream& operator<<(std::ostream& stream, const Constraint& c) {
	return stream << c.segNum << '\t'
				  << c.genome1 << '\t' << c.interval1 << '\t'
				  << c.genome2 << '\t' << c.interval2 << '\n';
}

void readAnchors(std::istream& stream, AnchorMap& anchors) {
	Anchor a;
	while (stream >> a) {
		anchors[a.id] = new Anchor(a);
	}
}

void makeClique(std::vector<std::string>& anchorNames,
				GenomeAnchorMaps& anchors,
				Clique& clique) {
	for (size_t i = 0; i < clique.size(); ++i) {
		std::string name = anchorNames[i];
		if (name == "NA") {
			clique[i] = NULL;
		} else {
			clique[i] = anchors[i][name];
		}
	}
}

bool isEmptyClique(const Clique& clique) {
	for (size_t i = 0; i < clique.size(); ++i) {
		if (clique[i] != NULL) { return false; }
	}
	return true;
}

void readRuns(std::istream& stream,
			  GenomeAnchorMaps& anchors,
			  RunMap& runs) {
	line::InputStream lineStream(stream);
	std::string line;
	std::vector<std::string> anchorNames;
	Clique clique(anchors.size());
	Run run;
	size_t num = 1;
	while (lineStream >> line) {
		anchorNames.clear();
		util::string::split(line, std::back_inserter(anchorNames));
		makeClique(anchorNames, anchors, clique);
		if (isEmptyClique(clique)) {
			runs[num] = run;
			++num;
			run.clear();
		} else {
			run.push_back(clique);
		}
	}
}

std::string underscoresToSpaces(const std::string& s) {
	std::string str = s;
	std::replace(str.begin(), str.end(), '_', ' ');
	return str;
}

genome::BasicInterval relativeInterval(const genome::Interval& insideInt,
									   const genome::Interval& enclosingInt) {
	genome::BasicInterval relativeInt;
	relativeInt.setChrom(insideInt.getChrom());
	if (enclosingInt.getStrand() == '+') {
		relativeInt.setStart(insideInt.getStart() - enclosingInt.getStart());
		relativeInt.setEnd(insideInt.getEnd() - enclosingInt.getStart());
		relativeInt.setStrand(insideInt.getStrand());
	} else {
		relativeInt.setStart(enclosingInt.getEnd() - insideInt.getEnd());
		relativeInt.setEnd(enclosingInt.getEnd() - insideInt.getStart());
		relativeInt.setStrand(insideInt.getStrand().opposite());
	}
	return relativeInt;
}

SimpleConstraint makeSimpleConstraint(const Constraint& c,
									  const genome::Interval& interval1,
									  const genome::Interval& interval2) {
	genome::BasicInterval relativeInt1 = relativeInterval(c.interval1, interval1);
	genome::BasicInterval relativeInt2 = relativeInterval(c.interval2, interval2);
	
	SimpleConstraint sc;
	sc.genome1 = c.genome1;
	sc.start1 = relativeInt1.getStart();
	sc.end1 = relativeInt1.getEnd();
	sc.genome2 = c.genome2;
	sc.start2 = relativeInt2.getStart();
	sc.end2 = relativeInt2.getEnd();
	return sc;
}

fasta::Record makeFASTARecord(SDB::DB& db,
							  const std::string& genome,
							  genome::Interval* interval) {
	return fasta::Record(genome, db.getSeq(*interval));
}

void writeClique(std::ostream& stream, Clique& c, homologymap::Segment& seg) {
	for (size_t i = 0; i < c.size(); ++i) {
		if (i > 0) { stream << '\t'; }
		if (c[i] == NULL) {
			stream << "-1\t-1";
		} else {
			genome::BasicInterval span = relativeInterval(c[i]->interval,
														  *seg.intervals[i]);
			stream << span.getStart() << '\t' << span.getEnd();
		}
	}
	stream << '\n';
}

void makeSegmentFiles(homologymap::Segment& seg,
					  phylogenetic::Tree* tree,
					  std::vector<std::string>& genomes,
					  unordered_map<std::string, size_t>& genomeNums,
					  std::vector<SDB::DB>& dbs,
					  std::vector<Constraint>& constraints,
					  Run& run,
					  fasta::OutputStream& softmaskedFile,
					  fasta::OutputStream& hardmaskedFile,
					  std::ostream& treeFile,
					  std::ostream& constraintsFile,
					  std::ostream& runFile) {
	unordered_set<std::string> includedGenomes;
	std::vector<fasta::Record> recs;

	for (size_t i = 0; i < seg.intervals.size(); ++i) {
		if (seg.hasGenome(i)) {
			includedGenomes.insert(underscoresToSpaces(genomes[i]));
			recs.push_back(makeFASTARecord(dbs[i], genomes[i], seg.intervals[i]));
		}
	}

	// Write sequence files
	for (size_t i = 0; i < recs.size(); ++i) {
		softmaskedFile << recs[i];
	}
	for (size_t i = 0; i < recs.size(); ++i) {
		bio::alphabet::Nucleotide::hardMaskInPlace(recs[i].sequence);
		hardmaskedFile << recs[i];
	}

	// Write tree file
	phylogenetic::Tree* subtree = tree->getSubtree(includedGenomes);
	newick::writeTree(treeFile, subtree);
	delete subtree;

	// Write constraints
	for (size_t i = 0; i < constraints.size(); ++i) {
		Constraint& c = constraints[i];
		constraintsFile <<
			makeSimpleConstraint(c,
								 *seg.intervals[genomeNums[c.genome1]],
								 *seg.intervals[genomeNums[c.genome2]]);
	}

	// Write run
	for (size_t i = 0; i < run.size(); ++i) {
		Clique& c = run[i];
		writeClique(runFile, c, seg);
	}
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Values to be (optionally) specified on the command line
	std::string mapDirname;
	std::string alignDirname;
	std::string mapFilename = "map";
	std::string treeFilename = "treefile";
	std::string softmaskedName = "seqs.fasta";
	std::string hardmaskedName = "seqs.fasta.masked";
	std::string treeName = "treefile";
	std::string constraintsName = "cons";
	std::string runsName = "runs";

	// Parse options
	util::options::Parser parser("",
								 "Create data files for a whole genome "
								 "alignment given an orthology map");
	parser.addStoreOpt(0, "map",
					   "Orthology map filename in map directory",
					   mapFilename, "FILE");
	parser.addStoreOpt(0, "tree",
					   "Tree filename in map directory",
					   treeFilename, "FILE");
	parser.addStoreOpt(0, "softmaskedName",
					   "Filename for softmasked sequence files",
					   softmaskedName);
	parser.addStoreOpt(0, "hardmaskedName",
					   "Filename for hardmasked sequence files",
					   hardmaskedName);
	parser.addStoreOpt(0, "treeName",
					   "Filename for tree files",
					   treeName);
	parser.addStoreOpt(0, "constraintsName",
					   "Filename for constraint files",
					   constraintsName);
	parser.addStoreOpt(0, "runsName",
					   "Filename for runs files",
					   runsName);
	parser.addStoreArg("mapDir",
					   "Directory containing map-related files",
					   mapDirname);
	parser.addStoreArg("alignDir",
					   "Diretory for alignment input files",
					   alignDirname);
	
	parser.parse(argv, argv + argc);

	try {
		Path mapDir(mapDirname);
		Path alignDir(alignDirname);
		
		// Read in genome names
		std::cerr << "Reading genome names...\n";
		std::vector<std::string> genomes;
		unordered_map<std::string, size_t> genomeNums;
		InputFileStream genomesFile(mapDir / "genomes");
		std::string genomesString = readStream(genomesFile);
		util::string::split(genomesString, std::back_inserter(genomes));
		for (size_t i = 0; i < genomes.size(); ++i) {
			genomeNums[genomes[i]] = i;
		}

		// Read tree file
		std::cerr << "Reading tree...\n";
		InputFileStream treeFile(mapDir / treeFilename);
		phylogenetic::Tree* tree = newick::readTree(treeFile);
		if (tree == NULL) {
			throw std::runtime_error("Invalid tree file: " + treeFilename);
		}

		// Create SDB objects for each genome
		std::cerr << "Opening SDB files...\n";
		std::vector<bio::SDB::DB> genomeDBs(genomes.size());
		for (size_t i = 0; i < genomes.size(); ++i) {
			genomeDBs[i].open(mapDir / (genomes[i] + ".sdb"));
			genomeDBs[i].readIndex();
		}

		// Read homology map
		std::cerr << "Reading map...\n";
		homologymap::Map map;
		InputFileStream mapFile(mapDir / mapFilename);
		map.read(mapFile);

		// Read anchors
		std::cerr << "Reading anchors...\n";
		GenomeAnchorMaps anchors(genomes.size());
		for (size_t i = 0; i < genomes.size(); ++i) {
			InputFileStream anchorsFile(mapDir / (genomes[i] + ".anchors"));
			readAnchors(anchorsFile, anchors[i]);
		}

		// Read runs
		std::cerr << "Reading runs...\n";
		RunMap runs;
		InputFileStream runsFile(mapDir / "runs");
		readRuns(runsFile, anchors, runs);

		// Read constraints
		std::cerr << "Reading constraints...\n";
		unordered_map<size_t, std::vector<Constraint> > constraints;
		InputFileStream constraintsFile(mapDir / "constraints");
		Constraint c;
		while (constraintsFile >> c) {
			constraints[c.segNum].push_back(c);
		}
		
		// Check validity of output path
		if (not alignDir.exists() or not alignDir.isDirectory()) {
			throw std::runtime_error("Invalid alignment directory: " +
									 alignDirname);
		}

		// Copy map and genomes file to alignment dir
		Path(mapDir / mapFilename).copyTo(alignDir / "map");
		Path(mapDir / "genomes").copyTo(alignDir / "genomes");

		// Create files for each segment
		for (size_t i = 0; i < map.getNumSegments(); ++i) {
			homologymap::Segment* seg = map.getSegment(i);

			// Create directory for this segment
			Path segDir = alignDir / toString(seg->num);
			if (not segDir.exists()) {
				segDir.createDirectory();
			}

			// Construct output files
			OutputFileStream softmaskedFile(segDir / softmaskedName);
			fasta::OutputStream softmaskedStream(softmaskedFile);
			OutputFileStream hardmaskedFile(segDir / hardmaskedName);
			fasta::OutputStream hardmaskedStream(hardmaskedFile);
			OutputFileStream treeFile(segDir / treeName);
			OutputFileStream constraintsFile(segDir/ constraintsName);
			OutputFileStream runsFile(segDir / runsName);

			// Make files
			std::cerr << seg->num << '\n';
			makeSegmentFiles(*seg,
							 tree,
							 genomes,
							 genomeNums,
							 genomeDBs,
							 constraints[seg->num],
							 runs[seg->num],
							 softmaskedStream,
							 hardmaskedStream,
							 treeFile,
							 constraintsFile,
							 runsFile);
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
