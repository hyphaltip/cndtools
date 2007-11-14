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

#include <cerrno>

#include "util/io.hh"
#include "filesystem/Path.hh"
#include "filesystem/InputFileStream.hh"
#include "filesystem/OutputFileStream.hh"

namespace filesystem {

	const std::string Path::SEPARATOR = "/";

	Path::Path(const char* s) : s(s) {}
	Path::Path(const std::string& s) : s(s) {}

	Path& Path::operator/=(const Path& p) {
		s += SEPARATOR;
		s += p.s;
		return *this;
	}
		
	Path Path::operator/(const Path& p) const {
		return Path(*this) /= p;
	}
	
	std::string Path::toString() const {
		return s;
	}

	Path Path::leaf() const {
		if (s.empty() or s == SEPARATOR) {
			return s;
		}
		
		std::string::size_type sepPos = s.rfind(SEPARATOR);
		if (sepPos == std::string::npos) {
			return s;
		} else {
			return s.substr(sepPos + 1);
		}
	}
	
	bool Path::exists() const {
		try {
			getStat();
		} catch (const Error& e) {
			if (e.getErrorNum() == ENOENT or e.getErrorNum() == ENOTDIR) {
				return false;
			} else {
				throw;
			}
		}
		return true;
	}
	
	bool Path::isDirectory() const {
		return S_ISDIR(getStat().st_mode);
	}

	bool Path::isRegularFile() const {
		return S_ISREG(getStat().st_mode);
	}
	
	bool Path::isSymbolicLink() const {
		return S_ISLNK(getStat().st_mode);
	}
	
	struct stat Path::getStat() const {
		struct stat path_stat;
		if (::stat(s.c_str(), &path_stat) != 0) {
			throw Error(s, errno);
		}
		return path_stat;
	}

    void Path::copyTo(const Path& target) const {
		struct stat source_stat(getStat());

		InputFileStream in(*this);
		
		// TODO: Use source permissions to open out file		
		OutputFileStream out(target);

		util::io::copy_stream(in, out);
	}
	
	void Path::createDirectory() const {
		if (::mkdir(s.c_str(), S_IRWXU|S_IRWXG|S_IRWXO) != 0) {
			throw Error(s, errno);
		}
	}

	void Path::openForInput(std::ifstream& strm) const {
		strm.open(s.c_str());
		if (not strm) {
			throw std::runtime_error("Could not open input file: " + s);
		}
	}
	
	void Path::openForOutput(std::ofstream& strm) const {
		strm.open(s.c_str());
		if (not strm) {
			throw std::runtime_error("Could not open output file: " + s);
		}
	}
	
	Path::Error::Error(const std::string& path, int errnum)
		: std::runtime_error(std::string(std::strerror(errnum)) + ": " + path),
		  path(path),
		  errnum(errnum)
	{}

	Path::Error::~Error() throw () {}
	
	std::string Path::Error::getPath() const {
		return path;
	}

	int Path::Error::getErrorNum() const {
		return errnum;
	}
	
}
