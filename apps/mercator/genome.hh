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

#ifndef __GENOME_HH__
#define __GENOME_HH__

#include "util/stl.hh"
using util::stl::hash_map;

#include "types.hh"

class Genome {
private:

	static vector<Genome*> genomes;
	static vector<Edge*> edges;
	static hash_map<string, Genome*> genomeMap;
	
	string name;
	size_t num;
	bool draft;

	vector<Chromosome*> chroms;

	hash_map<string, Chromosome*> chromMap;
	hash_map<string, Anchor*> anchorMap;

	typedef vector<Chromosome*>::iterator ChromIt;
	
	void loadChromosomeFile(const Path& dataDir);
	void loadAnchorFile(const Path& dataDir);
	void loadHitFile(const Path& dataDir, Genome* other,
					 const double maxEValue);
	void writeAGPFile(const Path dataDir);
	void writeAnchorFile(const Path dataDir);
	
	static void writeGenomePixelizerSetup(ostream& strm,
										  const GenomicDist sizeUnits);
	static void writeGenomePixelizerCoords(ostream& strm,
										  const GenomicDist sizeUnits);
	static void writeGenomePixelizerMatrix(ostream& strm);
	
 public:
	
	static void addGenome(string name, bool draft = false);
	static size_t getNumGenomes() { return genomes.size(); }
	static Genome* getGenome(const size_t g) { return genomes[g]; }
	static const vector<Genome*>& getGenomes() { return genomes; }
	static void loadFiles(const Path& dataDir,
						  const double maxEValue,
						  const std::string& phitFilename);
	static void loadPhits(const Path& phitFilename);
	static void initEdges(const float prunePct);
	static void writeGenomePixelizerFiles(const Path& dataDir);
	static string getGenomeNamesString();
	static string getGenomePixelizerSetupFilename();	
	static string getGenomePixelizerCoordsFilename();
	static string getGenomePixelizerMatrixFilename();
	static void writeAGPFiles(const Path& dataDir);
	static void writeAnchorFiles(const Path& dataDir);
	static void writeRunPermFiles(const Path& dataDir);

	static void assembleDraftGenomes();

	static void writeActiveHits(ostream& strm,
								size_t g1,
								size_t g2);

	const string& getName() const { return name; }
	size_t getNum() const { return num; }
	bool isDraft() const { return draft; }
	
	GenomicDist getLength() const;
	GenomicDist getCoverage() const;
	float getPctCoverage() const;

	void removeChrom(const Chromosome* chrom);
	size_t getNumChroms() const { return chroms.size(); }
	Chromosome* getChrom(const size_t c) const { return chroms[c]; }

	// Order chromosomes in this genome according to size and assign
	// global genomic coordinate offsets to each chromosome
	void makeGenomicCoords();

	// Make assembled chromosomes
	void assemble();

	// Unflip all anchors on this genome
	void unflipAnchors();
	
	size_t getNumAnchors() const;
	size_t getNumAnchorsInCliques() const;
	size_t getNumAnchorsInRuns() const;
	size_t getNumCliques() const;
	size_t getNumRuns() const;
	size_t getNumAnchorsRepetitive() const;
	
	void printCoverage(ostream& strm) const;	

	void markRepeats(const size_t repeatNum,
					 const float repeatPct);
	//	void filterRepeats();

	void writeRunPerm(std::ostream& strm) const;
	
};

#endif // __GENOME_HH__
