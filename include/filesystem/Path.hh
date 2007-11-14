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

#ifndef __FILESYSTEM_PATH_HH__
#define __FILESYSTEM_PATH_HH__

#include <string>
#include <stdexcept>
#include <fstream>
#include <sys/stat.h>

namespace filesystem {

    class Path {
	public:
		Path(const char* s = "");
		Path(const std::string& s);

		Path& operator/=(const Path& p);
		Path operator/(const Path& p) const;
		
		std::string toString() const;

		Path leaf() const;

		bool exists() const;
		bool isDirectory() const;
		bool isRegularFile() const;
		bool isSymbolicLink() const;
		
		void createDirectory() const;

		void copyTo(const Path& target) const;

		void openForInput(std::ifstream& strm) const;
		void openForOutput(std::ofstream& strm) const;
		
		class Error : public std::runtime_error {
		public:
			Error(const std::string& path, int errnum);
			~Error() throw ();
			std::string getPath() const;
			int getErrorNum() const;
		private:
			std::string path;
			int errnum;
		};
				
	private:
		struct stat getStat() const; 
		
		std::string s;
		static const std::string SEPARATOR;
    };

}

#endif // __FILESYSTEM_PATH_HH__
