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

#include "util/options.hh"
#include "filesystem/Path.hh"

namespace util {

	namespace options {

		std::string Argument::getHelp(size_t argColWidth,
									  size_t lineWidth) const {
			std::string help(2, ' ');
			help.append(argName);
		
			if (help.length() < argColWidth) {
				help = util::string::ljust(help, argColWidth);
			} else {
				help.append(4, ' ');
			}

			if (description.empty()) {
				return util::string::fill(help,
										argColWidth,
										std::string(2, ' '),
										std::string(2, ' '));
			} else {
				return util::string::fill(description,
										lineWidth,
										help,
										std::string(argColWidth, ' '));
			}
		}

		std::string Argument::getUsage() const {
			std::string usage;
			if (minNumArgs() == maxNumArgs()) {
				for (size_t i = minNumArgs(); i <= maxNumArgs(); ++i) {
					if (i != minNumArgs()) {
						usage += " ";
					}
					usage += argName;
				}
			} else if (minNumArgs() == 0 &&
					   maxNumArgs() == std::numeric_limits<size_t>::max()) {
				usage = "[" + argName + "]*";
			} else {
				if (minNumArgs() == 0) {
					usage = "[" + argName + "]";
				} else {
					usage = argName;
				}
				usage += "{";
				if (minNumArgs() != 0) {
					usage += util::string::toString(minNumArgs());
				}
				usage += ",";
				if (maxNumArgs() != std::numeric_limits<size_t>::max()) {
					usage += util::string::toString(maxNumArgs());
				}
				usage += "}";
			}
			return usage;
		}		
	
		std::string Option::getName() const {
			return shortName ?
				"-" + util::string::toString(shortName) :
				"--" + longName;
		}

		std::string Option::getHelp(size_t argColWidth,
									size_t lineWidth) const {
			std::string help(2, ' ');

			if (shortName) {
				help.append("-");
				help.append(1, shortName);
			} else {
				help.append(4, ' ');
			}

			if (!longName.empty()) {
				if (shortName) {
					help.append(", ");
				}
				help.append("--");
				help.append(longName);
			}

			if (takesOptionalArg() || numRequiredArgs() > 0) {
				if (longName.empty()) {
					help.append(argName);
				} else {
					help.append("=");
					help.append(argName);
				}
			}

			if (help.length() < argColWidth) {
				help = util::string::ljust(help, argColWidth);
			} else {
				help.append(4, ' ');
			}
			
			return util::string::fill(getDescription(),
									lineWidth,
									help,
									std::string(argColWidth, ' '));
		}



		Parser::Parser(const std::string& usage,
					   const std::string& progDescription)
			: usage(usage), progDescription(progDescription), progName(),
			  options(), arguments() {
			addOpt(new HelpOption('?', "help"));
		}
	
		void Parser::printHelp(std::ostream& strm) const {
			const size_t lineWidth = (getenv("COLUMNS") == NULL ?
									  80 : atoi(getenv("COLUMNS")));
			const size_t argColWidth = 30;
		
			printUsage(strm);
		
			if (!progDescription.empty()) {
				strm << util::string::fill(progDescription, lineWidth) << '\n';
			}

			if (!arguments.empty()) {
				strm << '\n' << "Arguments: " << '\n';
				for (std::vector<Argument*>::const_iterator it = arguments.begin();
					 it != arguments.end(); ++it) {
					if (!(*it)->description.empty()) {
						strm << (*it)->getHelp(argColWidth, lineWidth) << '\n';
					}
				}
			}
		
			if (!options.empty()) {
				strm << '\n' << "Options: " << '\n';
				for (std::vector<Option*>::const_iterator it = options.begin();
					 it != options.end(); ++it) {
					strm << (*it)->getHelp(argColWidth, lineWidth) << '\n';
				}
			}
		}

		void Parser::printUsage(std::ostream& strm) const {
			strm << "Usage: " << progName << " [options]";
			for (std::vector<Argument*>::const_iterator it = arguments.begin();
				 it != arguments.end(); ++it) {
				if (!(*it)->argName.empty()) {
					strm << ' ' << (*it)->getUsage();
				}
			}
			strm << ' ' << usage << '\n';
		}

		Option* Parser::findOption(char shortName) const {
			for (std::vector<Option*>::const_iterator it = options.begin();
				 it != options.end(); ++it) {
				if ((*it)->shortName == shortName) {
					return *it;
				}
			}
			throw OptionError("Invalid option -" + util::string::toString(shortName));
		}

		Option* Parser::findOption(const std::string& longName) const {
			Option* opt = NULL;
			for (std::vector<Option*>::const_iterator it = options.begin();
				 it != options.end(); ++it) {
				if ((*it)->longName == longName) {
					return *it;
				} else if (util::string::startsWith((*it)->longName, longName)) {
					if (opt == NULL) {
						opt = *it;
					} else {
						throw OptionError("Ambiguous option --" + longName);
					}
				}
			}
			if (opt == NULL) {
				throw OptionError("Invalid option --" + longName);
			}
			return opt;
		}

		bool Parser::isEndOfOptionsFlag(const std::string& str) {
			return str == "--";
		}
	
		bool Parser::isShortOption(const std::string& str) {
			return str.length() > 1 && str[0] == '-' && str[1] != '-';
		}
	
		bool Parser::isLongOption(const std::string& str) {
			return str.length() > 2 && str[0] == '-' && str[1] == '-';
		}
	
		void Parser::parse(const char** begin,
						   const char** end,
						   bool handleErrors,
						   std::ostream& errorStream) {

			// Extract program name
			assert(begin != end);
			filesystem::Path progPath(*begin);
			progName = progPath.leaf().toString();
			++begin;

			try {
				std::vector<Argument*>::iterator argIt = arguments.begin();
				size_t argCount = 0;
				bool processOptions = true;

				// Parse options and arguments
				while (begin != end) {
					std::string arg = *begin;
					++begin;

					if (processOptions && isEndOfOptionsFlag(arg)) {
						processOptions = false;
					} else if (processOptions && isShortOption(arg)) {
						parseShortOption(arg.substr(1), begin, end);
					} else if (processOptions && isLongOption(arg)) {
						parseLongOption(arg.substr(2), begin, end);
					} else if (argIt != arguments.end()) {
						(*argIt)->process(arg);
						++argCount;
						if ((*argIt)->maxNumArgs() == argCount) {
							argCount = 0;
							++argIt;
						}
					} else {
						throw ArgError("Wrong number of arguments given");
					}
				}

				// Check that the correct number of arguments was given
				if (argIt != arguments.end() &&
					argCount < (*argIt)->minNumArgs()) {
					throw ArgError("Wrong number of arguments given");
				}
				
			} catch (const HelpException& e) {
				if (!handleErrors) {
					throw e;
				}
				printHelp(errorStream);
				exit(EXIT_SUCCESS);
			} catch (const OptionError& e) {
				if (!handleErrors) {
					throw e;
				}
				errorStream << "Error: " << e.what() << '\n';
				exit(EXIT_FAILURE);
			} catch (const ArgError& e) {
				if (!handleErrors) {
					throw e;
				}
				errorStream << "Error: " << e.what() << '\n';
				printUsage(errorStream);
				exit(EXIT_FAILURE);
			}

		}

		void Parser::parseLongOption(const std::string& arg,
									 const char**& begin,
									 const char**& end) const {
			std::vector<std::string> optArgs;
			std::string::const_iterator equalSignPos =
				std::find(arg.begin() + 1, arg.end(), '=');
			std::string name(arg.begin(), equalSignPos);
			Option* opt = findOption(name);
			if (equalSignPos != arg.end()) {
				optArgs.push_back(std::string(equalSignPos + 1, arg.end()));
			}
			collectOptionArgs(opt, optArgs, begin, end);
			opt->process(optArgs);
		}

		void Parser::collectOptionArgs(Option* opt,
									   std::vector<std::string>& optArgs,
									   const char**& begin,
									   const char**& end) const {
			size_t argsAdded = 0;
			if (!opt->takesOptionalArg() && optArgs.size() == 1) {
				argsAdded = 1;
			}
			while (argsAdded < opt->numRequiredArgs() && begin != end) {
				optArgs.push_back(*begin);
				++begin;
				++argsAdded;
			}
			if (argsAdded != opt->numRequiredArgs()) {
				throw OptionError("Wrong number of arguments given for option "
								  + opt->getName());
			}
		}			
		
		void Parser::parseShortOption(const std::string& arg,
									  const char**& begin,
									  const char**& end) const {
			std::vector<std::string> optArgs;
			for (std::string::const_iterator c = arg.begin();
				 c != arg.end(); ++c) {
				Option* opt = findOption(*c);
				if (!opt->takesOptionalArg() &&
					opt->numRequiredArgs() == 0) {
					opt->process(optArgs);
				} else {
					++c;
					if (c != arg.end()) {
						optArgs.push_back(std::string(c, arg.end()));
					}
					collectOptionArgs(opt, optArgs, begin, end);
					opt->process(optArgs);
					break;
				}
			}
			
		}

	};	
};
