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

#ifndef __BIO_FORMATS_AGP_CONSTANTS_HH__
#define __BIO_FORMATS_AGP_CONSTANTS_HH__

#include <string>

namespace bio { namespace formats { namespace agp {

    class Constants {
	public:
		// Component types
		static const char ACTIVE_FINISHING;
		static const char DRAFT_HTG;
		static const char FINISHED_HTG;
		static const char WHOLE_GENOME_FINISHING;
		static const char PRE_DRAFT;
		static const char GAP;
		static const char OTHER;
		static const char WGS_CONTIG;

		// Gap types
		static const std::string FRAGMENT;
		static const std::string CONTIG;
		static const std::string SPLIT_FINISHED;
		static const std::string CLONE;
		static const std::string CENTROMERE;
		static const std::string SHORT_ARM;
		static const std::string HETEROCHROMATIN;
		static const std::string TELOMERE;

	protected:
		// Linkage strings
		static const std::string LINKAGE_NO;
		static const std::string LINKAGE_YES;
    };

} } }

#endif // __BIO_FORMATS_AGP_CONSTANTS_HH__
