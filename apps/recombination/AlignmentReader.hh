#include <iosfwd>
#include <string>

#include "bio/alignment/BasicNamedMultipleAlignment.hh"

class AlignmentReader {
public:
	AlignmentReader(std::istream& stream,
					const std::string& recombinant_name,
					size_t recombinant_num);

	const bio::alignment::MultipleAlignment& getAlignment() const;
	
	size_t getRecombinantNum() const;

private:
	bio::alignment::BasicNamedMultipleAlignment alignment;
	size_t recombinant_num;
};

inline
const bio::alignment::MultipleAlignment& AlignmentReader::getAlignment() const {
	return alignment;
}

inline
size_t AlignmentReader::getRecombinantNum() const { return recombinant_num; }
