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

#include "util/stl.hh"
#include "util/string.hh"
#include "util/io/line/InputStream.hh"

#include "genome.hh"
#include "anchor.hh"
#include "assembled.hh"
#include "clique.hh"
#include "chromosome.hh"
#include "edge.hh"
#include "run.hh"

vector<Genome*> Genome::genomes = vector<Genome*>();
vector<Edge*> Genome::edges = vector<Edge*>();
hash_map<string, Genome*> Genome::genomeMap = hash_map<string, Genome*>();

void Genome::unflipAnchors() {
	for_each(chroms.begin(), chroms.end(), mem_fun(&Chromosome::unflipAnchors));
}

void Genome::writeActiveHits(ostream& strm,
							 size_t g1,
							 size_t g2) {
	for (vector<Edge*>::iterator it = edges.begin();
		 it != edges.end(); ++it) {
		Edge* e = *it;
		if (e->hasGenomes(g1, g2) && e->isActive()) { 
			strm << e->getAnchor1()->getName() << '\t'
				 << e->getAnchor2()->getName() << '\t'
				 << e->getScore() << '\t'
				 << (e->getAnchor1()->isMarked() || e->getAnchor2()->isMarked() ?
					 '*' : ' ')
				 << '\n';
		}
	}
}

void Genome::addGenome(string name, bool draft) {
	Genome* genome = new Genome();
	genome->name = name;
	genome->num = genomes.size();
	genome->draft = draft;
	genomes.push_back(genome);
	genomeMap[name] = genome;
}

void Genome::initEdges(const float prunePct) {
	// 	std::for_each(edges.begin(), edges.end(), std::mem_fun(&Edge::addEdge));
	
	vector<Edge*>::iterator ePos;
	for (ePos = edges.begin(); ePos != edges.end(); ++ePos) {
		if (!(*ePos)->shouldBePruned(prunePct)) {
			(*ePos)->addEdge();
		}
	}
	
	// 	std::for_each(genomes.begin(), genomes.end(),
	// 				  std::bind2nd(std::mem_fun(&Genome::markRepeats), maxHits));
	// 	std::for_each(genomes.begin(), genomes.end(),
	// 				  std::mem_fun(&Genome::filterRepeats));
}

void Genome::markRepeats(const size_t repeatNum,
						 const float repeatPct) {
	vector<Chromosome*>::iterator it;
	for (it = chroms.begin(); it != chroms.end(); ++it) {
		for (Anchor* a = (*it)->getFirstAnchor();
			 a != NULL; a = a->nextAnchor()) {
			a->setMarked(a->isRepetitive(repeatNum, repeatPct));
		}
	}
}

// void Genome::filterRepeats() {
// 	vector<Chromosome*>::iterator it;
// 	for (it = chroms.begin(); it != chroms.end(); ++it) {
// 		for (Anchor* a = (*it)->getFirstAnchor();
// 			 a != NULL; a = a->nextAnchor()) {
// 			if (a->isMarked()) {
// 				a->removeAllEdges();
// 			}
// 		}
// 	}
// }

size_t Genome::getNumAnchors() const {
	return util::stl::sum(chroms.begin(), chroms.end(),
						  std::mem_fun(&Chromosome::getNumAnchors));
}

size_t Genome::getNumAnchorsInCliques() const {
	return util::stl::sum(chroms.begin(), chroms.end(),
						  std::mem_fun(&Chromosome::getNumAnchorsInCliques));
}

size_t Genome::getNumAnchorsInRuns() const {
	return util::stl::sum(chroms.begin(), chroms.end(),
						  std::mem_fun(&Chromosome::getNumAnchorsInRuns));
}

size_t Genome::getNumAnchorsRepetitive() const {
	return util::stl::sum(chroms.begin(), chroms.end(),
						  std::mem_fun(&Chromosome::getNumAnchorsRepetitive));
}

size_t Genome::getNumCliques() const {
	return util::stl::sum(chroms.begin(), chroms.end(),
						  std::mem_fun(&Chromosome::getNumCliques));
}

size_t Genome::getNumRuns() const {
	return util::stl::sum(chroms.begin(), chroms.end(),
						  std::mem_fun(&Chromosome::getNumRuns));
}

GenomicDist Genome::getLength() const {
	return util::stl::sum(chroms.begin(), chroms.end(),
						  std::mem_fun(&Chromosome::getLength));
}

GenomicDist Genome::getCoverage() const {
	return util::stl::sum(chroms.begin(), chroms.end(),
						  std::mem_fun(&Chromosome::getCoverage));
}

float Genome::getPctCoverage() const {
	if (getLength() == 0) {
		return 0;
	} else {
		return static_cast<float>(getCoverage()) / getLength();
	}
}

void Genome::removeChrom(const Chromosome* chrom) {
	chroms.erase(remove(chroms.begin(), chroms.end(), chrom),
				 chroms.end());
}

void Genome::makeGenomicCoords() {
	// Sort chromosomes by decreasing length
	sort(chroms.rbegin(), chroms.rend(), ChromLengthSorter());
	
	// Calculate genomic coordinates for each chrom
	vector<Chromosome*>::iterator chromPos;
	GenomicDist currStart = 0;
	for (chromPos = chroms.begin(); chromPos != chroms.end(); ++chromPos) {
		(*chromPos)->setGenomeStart(currStart);
		currStart += (*chromPos)->getLength();
	}
}

void Genome::printCoverage(ostream& strm) const {
	vector<Chromosome*>::const_iterator cPos;
	for (cPos = chroms.begin(); cPos != chroms.end(); ++cPos) {
		strm << (*cPos)->getName() << '\t'
			 << (*cPos)->getLength() << '\t'
			 << (*cPos)->getCoverage() << '\t'
			 << (*cPos)->getPctCoverage()
			 << '\n';
	}
	
	strm << "TOTAL" << '\t'
		 << getLength() << '\t'
		 << getCoverage() << '\t'
		 << getPctCoverage()
		 << '\n';
}

void Genome::loadChromosomeFile(const Path& dataDir) {
	// Construct path to chromosome file and check for existence
	Path chromPath = dataDir / (name + ".chroms");
	if (not chromPath.exists()) {
		throw std::runtime_error("Chromosome file does not exist: " +
								 chromPath.toString());
	}

	// Read in chromosome file
	std::cerr << name;
	InputFileStream chromFile(chromPath);
	while (chromFile) {
		std::string name = "";
		GenomicDist length = 0;
		chromFile >> name >> length;
		if (not name.empty()) {
			Chromosome* chrom = new Chromosome(this, name, length);
			chroms.push_back(chrom);
			chromMap[name] = chrom;
		}
	}

	makeGenomicCoords();

	// Note how many chromosomes were read
	std::cerr << " " << getNumChroms() << " "
			  << (isDraft() ? "contigs" : "chromosomes") << '\n';
}

void Genome::loadAnchorFile(const Path& dataDir) {
	// Construct path to anchor file and check for existence
	Path infoPath = dataDir / (name + ".anchors");
	if (not infoPath.exists()) {
		throw std::runtime_error("Anchor file does not exist: " +
								 infoPath.toString());
	}
	
	// Read in info file
	std::cerr << name;
	InputFileStream infoFile(infoPath);
   	
	while (infoFile) {
		string aname =  "";
		string chrom = "";
		string strand = "";
		GenomicDist start = 0;
		GenomicDist end = 0;
		size_t isCoding = 0;

		infoFile >> aname >> chrom >> strand >> start >> end >> isCoding;
		
		if (not aname.empty()) {
			if (chromMap.find(chrom) == chromMap.end()) {
				cerr << "WARNING: anchor '" << aname << "'"
					 << " is in a chromosome '" << chrom << "'"
					 << " not listed in the chromosome file" << endl;
				continue;
			}

			Chromosome* c = chromMap[chrom];
			Anchor* a = new Anchor(aname, c, strand.at(0), start, end, isCoding);
			c->addAnchor(a);
			
			anchorMap[aname] = a;
		}
	}

	// Initialize (sort, etc.) anchors in each chromosome
	std::for_each(chroms.begin(), chroms.end(),
				  std::mem_fun(&Chromosome::initAnchors));

	// Note how many anchors are in this genome
	std::cerr << " " << getNumAnchors() << " anchors" << '\n';
}

void Genome::loadHitFile(const Path& dataDir,
						 Genome* other,
						 const double maxEValue) {
	// Construct two possible paths
	Genome* g1 = this;
	Genome* g2 = other;
	Path hitPath = dataDir / (g1->getName() + "-" + g2->getName() + ".hits");
	if (not hitPath.exists()) {
		std::swap(g1, g2);
		hitPath = dataDir / (g1->getName() + "-" + g2->getName() + ".hits");
		if (not hitPath.exists()) {
			throw std::runtime_error("Hit file for " + this->getName() +
									 " and " + other->getName() +
									 " could not be found");
		}
	}

	// Read in hit file
	std::cerr << (g1->getName() + "-" + g2->getName());
	InputFileStream hitFile(hitPath);
	util::io::line::InputStream hitStream(hitFile);
	
	size_t numEdgesBefore = edges.size();
	size_t numEdgesFiltered = 0;
	
	std::string name1;
	std::string name2;
	int score;
	double evalue;

	std::string line;
	while (hitStream >> line) {
		std::istringstream lineStream(line);
		lineStream >> name1 >> name2 >> score >> evalue;

		if (not lineStream) {
			throw std::runtime_error("Invalid hit: " + line);
		}
			
		Anchor* a1 = g1->anchorMap[name1];
		Anchor* a2 = g2->anchorMap[name2];

		if (a1 == NULL || a2 == NULL) {
			std::cerr << "WARNING: Found hit including an anchor not listed "
					  << "in the anchor files: "
					  << name1 << " " << name2 << '\n';
			continue;
		}
		
		if (evalue <= maxEValue) {
			Edge* e = new Edge(a1, a2, score);
			edges.push_back(e);
		} else {
			++numEdgesFiltered;
		}
	}
	std::cerr << " " << edges.size() - numEdgesBefore << " hits"
			  << " (" << numEdgesFiltered << " filtered)" << '\n';
}

void Genome::loadPhits(const Path& phitFilename) {
	InputFileStream phitFile(phitFilename);

	int segNum;
	std::string gname1;
	std::string gname2;
	std::string aname1;
	std::string aname2;
	
	while (phitFile >> segNum >> gname1 >> aname1 >> gname2 >> aname2) {
		Genome* g1 = genomeMap[gname1];
		Genome* g2 = genomeMap[gname2];
		if (g1 == NULL) {
			throw std::runtime_error("Unspecified genome in phitfile: " +
									 gname1);
		}
		if (g2 == NULL) {
			throw std::runtime_error("Unspecified genome in phitfile: " +
									 gname2);
		}
		
		Anchor* a1 = g1->anchorMap[aname1];
		Anchor* a2 = g2->anchorMap[aname2];

		if (a1 == NULL || a2 == NULL) {
			std::cerr << "WARNING: Found hit including an anchor not listed "
					  << "in the anchor files: "
					  << aname1 << " " << aname2 << '\n';
			continue;
		}
		
		Edge* e = new Edge(a1, a2, 1);
		edges.push_back(e);
	}
	
	std::cerr << " " << edges.size() << " pairwise hits loaded\n";
}
	

void Genome::loadFiles(const Path& dataDir,
					   const double maxEValue,
					   const std::string& phitFilename) {
	vector<Genome*>::iterator gPos;
	
	// First load all of the chromosome files
	cerr << "Loading chromosome files...\n";
	for (gPos = genomes.begin(); gPos != genomes.end(); ++gPos) {
		(*gPos)->loadChromosomeFile(dataDir);
	}
	
	// Load all of the anchor files
	cerr << "Loading anchor files...\n";
	for (gPos = genomes.begin(); gPos != genomes.end(); ++gPos) {
		(*gPos)->loadAnchorFile(dataDir);
	}	

	if (phitFilename != "") {
		cerr << "Loading pairwise hits...\n";
		loadPhits(phitFilename);
	} else {
		// Now load hit files
		cerr << "Loading hit files...\n";
		for (gPos = genomes.begin(); gPos != genomes.end(); ++gPos) {
			vector<Genome*>::iterator gPos2;
			for (gPos2 = gPos + 1; gPos2 != genomes.end(); ++gPos2) {
				(*gPos)->loadHitFile(dataDir, *gPos2, maxEValue);
			}
		}

		cerr << "Sorting edges...\n";
		sort(edges.rbegin(), edges.rend(), EdgeSorter());
	}

}

string Genome::getGenomeNamesString() {
	string genomeStr = genomes[0]->getName();
	for (size_t g = 1; g < getNumGenomes(); ++g) {
		genomeStr += "-" + genomes[g]->getName();
	}
	return genomeStr;
}

string Genome::getGenomePixelizerSetupFilename() {
	return getGenomeNamesString() + ".genpix.setup";
}

string Genome::getGenomePixelizerCoordsFilename() {
	return getGenomeNamesString() + ".genpix.coords";
}

string Genome::getGenomePixelizerMatrixFilename() {
	return getGenomeNamesString() + ".genpix.matrix";
}		

void Genome::writeGenomePixelizerSetup(ostream& strm,
									   const GenomicDist sizeUnits) {
	// Get the sizes of all the chromosomes
	vector<GenomicDist> chromSizes;
	for (vector<Genome*>::const_iterator gpos = genomes.begin();
		 gpos != genomes.end(); ++gpos) {
		for (vector<Chromosome*>::const_iterator cpos = (*gpos)->chroms.begin();
			 cpos != (*gpos)->chroms.end(); ++cpos) {
			chromSizes.push_back((*cpos)->getLength());
		}
	}

	// Write setup information
	
	strm << "1. name of file containing gene coordinates: "
		 << getGenomePixelizerCoordsFilename() << '\n'
		 << "2. name of the distance matrix file: "
		 << getGenomePixelizerMatrixFilename() << '\n'
		 << "3. number of chromosomes: "
		 << chromSizes.size() << '\n'
		 << "4. size of chromosomes:";

	// write out chromosome sizes
	for (vector<GenomicDist>::const_iterator spos = chromSizes.begin();
		 spos != chromSizes.end(); ++spos) {
		strm << ' ' << static_cast<float>(*spos) / sizeUnits;
	}
	strm << '\n';

	strm << "5. identity upper level: 100" << '\n'
		 << "6. identity lower level: 97" << '\n'
		 << "7. window size (pixels) X: 800" << '\n'
		 << "8. window size (pixels) Y: 600" << '\n'
		 << "9. html prefix: http://bio.math.berkeley.edu" << '\n'
		 << "10. Title: Project Name" << '\n'
		 << "11. Laboratory: Pachter Lab" << '\n'
		 << "########################################################" << '\n'
		 << "#####   for experienced users below this line   ########" << '\n'
		 << "12. W/C correction: A" << '\n'
		 << "13. horizontal size of gene: 12" << '\n'
		 << "14. vertical size of gene: 6" << '\n'
		 << "15. W/C coefficient: 1" << '\n'
		 << "16. W/C correction value: 8" << '\n'
		 << "17. chromosome thickness: 5" << '\n'
		 << "18. gene feature mode (standard [std] or extended [ext]): std"
		 << '\n';
}

void Genome::writeGenomePixelizerCoords(ostream& strm,
										const GenomicDist sizeUnits) {
	size_t chromNum = 1;
	for (size_t g = 0; g < genomes.size(); ++g) {
		for (size_t c = 0; c < genomes[g]->chroms.size(); ++c) {
			Chromosome* chrom = genomes[g]->chroms[c];
			for (Anchor* a = chrom->getFirstAnchor();
				 a != NULL; a = a->nextAnchor()) {
				strm << chromNum << '\t'
					 << a->getUniqueName() << '\t'
					 << static_cast<float>(a->getStart()) / sizeUnits << '\t'
					 << (a->getStrand() == '+' ? 'W' : 'C') << '\t'
					 << "orange"
					 << '\n';
			}
			chromNum += 1;
		}
	}
}

void Genome::writeGenomePixelizerMatrix(ostream& strm) {
	for (size_t i = 0; i < edges.size(); ++i) {
		Edge* e = edges[i];
		strm << e->getAnchor1()->getUniqueName() << '\t'
			 << e->getAnchor2()->getUniqueName() << '\t'
			 << e->getScore()
			 << '\n';
	}
}

void Genome::writeGenomePixelizerFiles(const Path& dataDir) {
	const GenomicDist sizeUnits = 1000000;
	OutputFileStream setupFile(dataDir / getGenomePixelizerSetupFilename());
	writeGenomePixelizerSetup(setupFile, sizeUnits);
	OutputFileStream coordsFile(dataDir / getGenomePixelizerCoordsFilename());
	writeGenomePixelizerCoords(coordsFile, sizeUnits);
	OutputFileStream matrixFile(dataDir / getGenomePixelizerMatrixFilename());
	writeGenomePixelizerMatrix(matrixFile);	
} 

void Genome::writeAnchorFiles(const Path& dataDir) {
	std::for_each(genomes.begin(), genomes.end(),
				  std::bind2nd(std::mem_fun(&Genome::writeAnchorFile), dataDir));
}

void Genome::writeAGPFiles(const Path& dataDir) {
	std::for_each(genomes.begin(), genomes.end(),
				  std::bind2nd(std::mem_fun(&Genome::writeAGPFile), dataDir));
}

void Genome::writeRunPermFiles(const Path& dataDir) {
	vector<Genome*>::iterator it;
	for (it = genomes.begin(); it != genomes.end(); ++it) {
		OutputFileStream permFile(dataDir / ((*it)->getName() + ".mgr"));
		(*it)->writeRunPerm(permFile);
	}
}

void Genome::writeAGPFile(const Path dataDir) {
	// Construct path to assembled chrom file
	OutputFileStream assemFile(dataDir / (name + ".agp"));

	size_t recNum = 1;
	for (ChromIt it = chroms.begin(); it != chroms.end(); ++it) {
		(*it)->writeAGPLines(assemFile, recNum);
	}
}

void Genome::writeAnchorFile(const Path dataDir) {
	// Construct path to assembled anchor file
	OutputFileStream assemFile(dataDir / (name + ".anchors"));

	for (ChromIt it = chroms.begin(); it != chroms.end(); ++it) {
		for (Anchor* a = (*it)->getFirstAnchor(); a != NULL; a = a->nextAnchor()) {
			a->writeAnchorLine(assemFile);
		}
	}
}

void Genome::assembleDraftGenomes() {
	vector<Genome*>::iterator it;
	for (it = genomes.begin(); it != genomes.end(); ++it) {
		if ((*it)->isDraft()) {
			(*it)->assemble();
			(*it)->makeGenomicCoords();
		}
	}
}

void Genome::assemble() {
	size_t aChromNum = 1;
	
	// Get all runs that include this genome
	vector<Run*> runs;
	getRuns(runs, num);

	vector<Run*>::iterator it;
	for (it = runs.begin(); it != runs.end(); ++it) {
		// Continue if this run has already been visited
		if ((*it)->wasVisited()) {
			continue;
		}
		
		// Traverse all the way to the left from this run
		Run* curr = (*it);
		for (Run* next = curr->leftRun(num); next != NULL;
			 curr = next, next = curr->leftRun(num)) {

			assert(!next->wasVisited());

			// Orient next run so that its right side matches up with
			// the left side of the current run
			if (curr->getLeftChrom(num) != next->getRightChrom(num) ||
				curr->isLeftForward(num) != next->isRightForward(num)) {
				next->flip();
				assert(curr->getLeftChrom(num) == next->getRightChrom(num) &&
					   curr->isLeftForward(num) == next->isRightForward(num));
			}
		}

		Assembled* assemChrom = new Assembled(this,
											  string("assembled") +
											  util::string::toString(aChromNum));
		++aChromNum;
		
		// Starting from the left-most run, add the chromosomes of the
		// cliques to a new assembled chromosome
		while (true) {
			curr->setVisited();
			
			vector<Clique*> cliques;
			curr->getCliques(cliques, num);

			vector<Clique*>::iterator clique_it;
			for (clique_it = cliques.begin(); clique_it != cliques.end();
				 ++clique_it) {
				Anchor* currAnchor = (*clique_it)->getAnchor(num);
				
				// Check if this anchor has already been moved onto
				// the assembled chromosome
				if (currAnchor->getChrom() == assemChrom) {
					assert(currAnchor->isForward());
					continue;
				} else {
					// Make sure anchors are all forward
					if (!currAnchor->isForward()) {
						currAnchor->getChrom()->reverse();
						assert(currAnchor->isForward());
					}
					// Add this anchor's chromosome to the assembled chromosome
					assemChrom->addChrom(currAnchor->getChrom());
				}
			}

			// Move to next run to the right
			Run* next = curr->rightRun(num);

			if (next == NULL) {
				break;
			}

			// Orient next run so that its left side matches up with
			// the right side of the current run
			if (curr->getRightChrom(num) != next->getLeftChrom(num) ||
				curr->isRightForward(num) != next->isLeftForward(num)) {
				next->flip();
				assert(curr->getRightChrom(num) == next->getLeftChrom(num) &&
					   curr->isRightForward(num) == next->isLeftForward(num));
			}

			curr = next;
		}

		chroms.push_back(assemChrom);
	}

	// reset visited flags
	for_each(runs.begin(), runs.end(),
			 bind2nd(mem_fun(&Run::setVisited), false));

	// remove chromosomes that were assembled into new chromosomes
	chroms.erase(remove_if(chroms.begin(), chroms.end(),
						   mem_fun(&Chromosome::isPartOfAssembled)),
				 chroms.end());
}

void Genome::writeRunPerm(std::ostream& strm) const {
	strm << ">" << name << '\n';
	for (size_t c = 0; c < chroms.size(); ++c) {
		chroms[c]->writeRunPerm(strm);
	}
}
