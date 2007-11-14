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
#include <vector>

#include "bio/sdb.hh"
#include "bio/genome/BasicInterval.hh"
#include "bio/formats/fasta.hh"
#include "bio/phylogenetic/Tree.hh"
#include "bio/formats/newick.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "util/stl.hh"
#include "filesystem.hh"
using namespace bio;
using namespace bio::genome;
using namespace bio::formats;
using namespace filesystem;
using util::string::toString;
using util::stl::hash_map;
using util::stl::hash_set;

struct Segment {
	size_t num;
	std::string genome;
	BasicInterval interval;
};

struct Edge {
	size_t num;
	size_t segNum1;
	size_t segNum2;
	Strand strand1;
	Strand strand2;
};

typedef hash_map<size_t, Segment> SegmentMap;
typedef std::vector<Edge> EdgeList;
typedef hash_map<std::string, SDB::DB> Genomes;

std::istream& operator>>(std::istream& stream, Segment& s) {
	size_t num;
	std::string genome;
	std::string chrom;
	Position start, end;
	stream >> num >> genome >> chrom >> start >> end;
	if (stream) {
		s.num = num;
		s.genome = genome;
		s.interval.setChrom(chrom);
		s.interval.setStart(start);
		s.interval.setEnd(end);
	}
	return stream;
}

std::istream& operator>>(std::istream& stream, Edge& e) {
	return stream >> e.num
				  >> e.segNum1 >> e.segNum2
				  >> e.strand1 >> e.strand2;
}

void readSegments(std::istream& stream, SegmentMap& segments) {
	Segment s;
	while (stream >> s) {
		segments[s.num] = s;
	}
}

void readEdges(std::istream& stream, EdgeList& edges) {
	Edge e;
	while (stream >> e) {
		edges.push_back(e);
	}
}

void openGenomes(Genomes& genomes, SegmentMap& segments, const Path& mapDir) {
	for (SegmentMap::iterator seg = segments.begin();
		 seg != segments.end(); ++seg) {
		std::string genome = seg->second.genome;
		if (genomes.find(genome) == genomes.end()) {
			SDB::DB& db = genomes[genome];
			db.open(mapDir / (genome + ".sdb"));
			db.readIndex();
		}
	}
}

std::string getSeq(Segment& s,
				   Strand strand,
				   SDB::DB& db) {
	BasicInterval interval(s.interval);
	interval.setStrand(strand);
	return db.getSeq(interval);
}

std::string underscoresToSpaces(const std::string& s) {
	std::string str = s;
	std::replace(str.begin(), str.end(), '_', ' ');
	return str;
}

void makeEdgeFiles(Edge& e,
				   SegmentMap& segments,
				   Genomes& genomes,
				   phylogenetic::Tree* tree,
				   std::ostream& softmaskedStream,
				   std::ostream& hardmaskedStream,
				   std::ostream& treeStream) {
	const std::string SEQ1_TITLE = "seq1";
	const std::string SEQ2_TITLE = "seq2";

	fasta::OutputStream softFASTAStream(softmaskedStream);
	fasta::OutputStream hardFASTAStream(hardmaskedStream);

	Segment& seg1 = segments[e.segNum1];
	Segment& seg2 = segments[e.segNum2];

	fasta::Record
		rec1(SEQ1_TITLE, getSeq(seg1, e.strand1, genomes[seg1.genome])),
		rec2(SEQ2_TITLE, getSeq(seg2, e.strand2, genomes[seg2.genome]));
	
	// Write sequence files
	softFASTAStream << rec1 << rec2;

	bio::alphabet::Nucleotide::hardMaskInPlace(rec1.sequence);
	bio::alphabet::Nucleotide::hardMaskInPlace(rec2.sequence);
	hardFASTAStream << rec1 << rec2;
	
	// Write tree file
	hash_set<std::string> includedGenomes;
	includedGenomes.insert(underscoresToSpaces(seg1.genome));
	includedGenomes.insert(underscoresToSpaces(seg2.genome));
	phylogenetic::Tree* subtree = tree->getSubtree(includedGenomes);
	if (subtree->getNumDescendants() != 2) {
		throw std::runtime_error("Tree does not contain genomes: " +
								 seg1.genome + " and/or " + seg2.genome);
	}
	subtree->getDescendant(0)->setLabel(SEQ1_TITLE);
	subtree->getDescendant(1)->setLabel(SEQ2_TITLE);
	newick::writeTree(treeStream, subtree);
	delete subtree;
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Values to be (optionally) specified on the command line
	std::string mapDirname = ".";
	std::string outDirname = ".";
	std::string edgeFilename = "edges";
	std::string segmentFilename = "segments";
	std::string treeFilename = "treefile";
	std::string genomesFilename = "genomes";
	std::string softmaskedName = "seqs.fasta";
	std::string hardmaskedName = "seqs.fasta.masked";
	std::string treeName = "treefile";

	// Parse options
	util::options::Parser parser("",
								 "Create data files for a whole genome "
								 "alignment given an orthology map");
	parser.addStoreOpt(0, "map-dir",
					   "Directory containing map-related files",
					   mapDirname, "DIR");
	parser.addStoreOpt(0, "out-dir",
					   "Output directory",
					   outDirname, "DIR");
	parser.addStoreOpt(0, "tree",
					   "Tree filename",
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
	parser.addStoreOpt(0, "edges",
					   "Edges filename",
					   edgeFilename);
	parser.addStoreOpt(0, "segments",
					   "Segments filename",
					   segmentFilename);
	parser.parse(argv, argv + argc);

	try {
		Path mapDir(mapDirname);
		Path outDir(outDirname);

		// Read segments
		std::cerr << "Reading segments...\n";
		InputFileStream segmentFile(mapDir / segmentFilename);
		SegmentMap segments;
		readSegments(segmentFile, segments);

		// Read edges
		std::cerr << "Reading edges...\n";
		InputFileStream edgeFile(mapDir / edgeFilename);
		EdgeList edges;
		readEdges(edgeFile, edges);

		// Get genome names
		std::cerr << "Opening genomes files...\n";
		Genomes genomes;
		openGenomes(genomes, segments, mapDir);

		// Read tree file
		std::cerr << "Reading tree...\n";
		InputFileStream treeFile(mapDir / treeFilename);
		phylogenetic::Tree* tree = newick::readTree(treeFile);
		if (tree == NULL) {
			throw std::runtime_error("Invalid tree file: " +
									 treeFilename);
		}
		
		// Check validity of output path
		if (not outDir.exists() or not outDir.isDirectory()) {
			throw std::runtime_error("Invalid output directory: " +
									 outDirname);
		}

		// Create files for each edge
		for (size_t i = 0; i < edges.size(); ++i) {
			Edge& e = edges[i];

			// Create directory for this edge
			Path edgeDir = outDir / toString(e.num);
			if (not edgeDir.exists()) {
				edgeDir.createDirectory();
			}
			
			// Construct output files
			OutputFileStream softmaskedFile(edgeDir / softmaskedName);
			OutputFileStream hardmaskedFile(edgeDir / hardmaskedName);
			OutputFileStream treeFile(edgeDir / treeName);
			
			makeEdgeFiles(e, segments, genomes, tree,
						  softmaskedFile, hardmaskedFile, treeFile);
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
