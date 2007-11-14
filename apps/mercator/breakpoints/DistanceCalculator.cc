#include <stdexcept>

#include "DistanceCalculator.hh"
#include "Genome.hh"

#include "bio/formats/newick.hh"

DistanceCalculator::DistanceCalculator(std::istream& strm) {
	readTree(strm);
}

bio::phylogenetic::Tree* DistanceCalculator::getTree() {
	return tree;
}
	
void DistanceCalculator::readTree(std::istream& strm) {
	// Read and parse tree file
	tree = bio::formats::newick::readTree(strm);
	if (tree == NULL) {
		throw std::runtime_error("Tree not in proper newick format");
	}
}

std::string DistanceCalculator::spacesToUnderscores(const std::string& s) {
	std::string str = s;
	std::replace(str.begin(), str.end(), ' ', '_');
	return str;
}

util::Matrix<double> DistanceCalculator::getDistances() {
	util::Matrix<double> distances(Genome::getNumGenomes(),
								   Genome::getNumGenomes());
	
	tree->setNumbers();
	typedef std::vector<const bio::phylogenetic::Tree*> TreeVect;
	TreeVect trees;
	tree->getTopologicalOrder(trees);
	util::Matrix<double> allDistances;
	tree->getPairwiseDistances(allDistances);

	for (size_t i = 0; i < allDistances.getNumRows(); ++i) {
		std::string label1 = spacesToUnderscores(trees[i]->getLabel());
		if (!trees[i]->isLeaf() || !Genome::getGenome(label1)) {
			continue;
		}
		for (size_t j = 0; j < allDistances.getNumCols(); ++j) {
			std::string label2 = spacesToUnderscores(trees[j]->getLabel());
			if (!trees[j]->isLeaf() || !Genome::getGenome(label2)) {
				continue;
			}
			size_t genome1 = Genome::getGenome(label1)->getNum();
			size_t genome2 = Genome::getGenome(label2)->getNum();
			distances(genome1, genome2) = allDistances(i, j);
		}
	}

	return distances;
}
