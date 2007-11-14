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

#include "bio/gff/GFFRecord.hh"
#include "bio/gff/GFFInputStream.hh"
#include "bio/homologymap/HomologyMapper.hh"
#include "util/string.hh"
#include "util/options.hh"
using util::string::toString;
using bio::gff::GFFRecord;
using bio::genome::BasicInterval;
using bio::homologymap::HomologyMapper;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool output_unmapped = false;
	std::string unmapped_attribute = "UNMAPPED";
	std::string seg_attr = "segment";
	std::string align_dir;	
	std::string source;
	std::string target;
	
	util::options::Parser parser("< gffInput", "");
	parser.addStoreOpt('s', "segattr",
					   "name of attribute to add to mapped features that "
					   "are broken up into multiple segments",
					   seg_attr, "ATTRIBUTE");
	parser.addStoreTrueOpt(0, "output-unmapped",
						   "Output unmapped records.  An extra attribute "
						   "(see --unmapped-attribute option) will be added "
						   "to unmapped records", output_unmapped);
	parser.addStoreOpt(0, "unmapped-attribute",
					   "Attribute to add to unmapped records",
					   unmapped_attribute);
	parser.addStoreArg("align_dir", "", align_dir);
	parser.addStoreArg("source_genome", "", source);
	parser.addStoreArg("target_genome", "", target);
	parser.parse(argv, argv + argc);
	
	try {
		HomologyMapper mapper(align_dir, source, target);

		bio::gff::GFFInputStream gffStream(std::cin);

		GFFRecord rec;
		std::vector<BasicInterval> mapped;
		while (gffStream >> rec) {
			mapped.clear();
			mapper.map(rec.getInterval(), mapped);

			if (mapped.empty()) {
				if (output_unmapped) {
					rec.addAttribute(unmapped_attribute);
					std::cout << rec;
					continue;
				} else {
					std::cerr << "Warning: Record not mapped: " << rec;
				}
			}
			
			if (not rec.hasAttribute(seg_attr)) {
				rec.addAttribute(seg_attr);
			}
			GFFRecord::Attribute& attr = rec.getAttribute(seg_attr);
			attr.values.push_back("0");
			std::string& seg_num = attr.values.back();
			for (size_t i = 0; i < mapped.size(); ++i) {
				rec.setInterval(mapped[i]);
				seg_num = toString(i + 1);
				std::cout << rec;
			}
		}
			
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
