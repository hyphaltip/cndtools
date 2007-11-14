#include <stdexcept>

#include "AlignmentReader.hh"
#include "bio/formats/fasta.hh"

AlignmentReader::AlignmentReader(std::istream& stream,
								 const std::string& recombinant_name,
								 size_t recombinant_num) {
	bio::formats::fasta::InputStream fastaStream(stream);
	fastaStream >> alignment;
	
	if (not recombinant_name.empty()) {
		int seq_num = alignment.getSeqNum(recombinant_name);
		if (seq_num == -1) {
			throw std::runtime_error("Invalid sequence name: " +
									 recombinant_name);
		}
		recombinant_num = seq_num;
	} else if (recombinant_num >= alignment.getNumSeqs()) {
		throw std::runtime_error("Invalid sequence number");
	}
	this->recombinant_num = recombinant_num;
}
