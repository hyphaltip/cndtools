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

#ifndef __BIO_FORMATS_MAF_CONSTANTS_HH__
#define __BIO_FORMATS_MAF_CONSTANTS_HH__

#include <string>

namespace bio { namespace formats { namespace maf {

    class Constants {
	protected:
		static const std::string VERSION;
		static const std::string ALIGNMENT_LINE_PREFIX;
		static const std::string SEQ_LINE_PREFIX;
		static const std::string HEADER_LINE_PREFIX;
		static const std::string METADETA_LINE_PREFIX;
		static const std::string COMMENT_LINE_PREFIX;
	};

} } }

#endif // __BIO_FORMATS_MAF_CONSTANTS_HH__
