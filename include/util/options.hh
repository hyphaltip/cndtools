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

#ifndef __OPTIONS_HH__
#define __OPTIONS_HH__

#include <string>
#include <stdexcept>
#include <vector>
#include <limits>
#include <iostream>

#include "boost/lexical_cast.hpp"

#include "util/string.hh"

namespace util {

	namespace options {

		struct HelpException : public std::exception {
		};

		struct ArgError : public std::runtime_error {
			ArgError(const std::string& what) : std::runtime_error(what) {}
		};

		struct OptionError : public std::runtime_error {
			OptionError(const std::string& what) : std::runtime_error(what) {}
		};
	
		class Argument {
		public:
			std::string argName;
			std::string description;

			Argument(const std::string& argName,
					 const std::string& description)
				: argName(argName), description(description) {
			}

			virtual ~Argument() {}
		
			// Invoke this option with arguments given by ARGS
			virtual void process(const std::string& val) { return; }

			virtual size_t minNumArgs() const { return 1; }
			virtual size_t maxNumArgs() const { return 1; }
		
			// Print the usage information for this argument to STRM
			virtual std::string getHelp(size_t argColWidth = 30,
										size_t lineWidth = 80) const;

			virtual std::string getUsage() const;
		};

		template<typename ArgType>
		class StoreArgument : public Argument {
		public:
			ArgType& arg;
		
			StoreArgument(const std::string& argName,
						  const std::string& description,
						  ArgType& arg)
				: Argument(argName, description), arg(arg) {
			}

			virtual void process(const std::string& val) {
				try {
					arg = boost::lexical_cast<ArgType, const std::string>(val);
				} catch (boost::bad_lexical_cast& e) {
					throw ArgError("Bad value for argument " + argName);
				}
			}
		};

		template<typename Container>
		class AppendArgument : public Argument {
		public:
			Container& c;
			size_t minArgs;
			size_t maxArgs;
		
			AppendArgument(const std::string& argName,
						   const std::string& description,
						   Container& c,
						   size_t minArgs = 0,
						   size_t maxArgs = std::numeric_limits<size_t>::max())
				: Argument(argName, description), c(c),
				  minArgs(minArgs), maxArgs(maxArgs) {
			}

			virtual size_t minNumArgs() const { return minArgs; }
			virtual size_t maxNumArgs() const { return maxArgs; }
		
			virtual void process(const std::string& val) {
				typedef typename Container::value_type ValueType;
				try {
					c.push_back(boost::lexical_cast<ValueType, const std::string>(val));
				} catch (boost::bad_lexical_cast& e) {
					throw ArgError("Bad value for argument " + argName);
				}
			}

		};
	
		class Option {
		public:
			char shortName;
			std::string longName;
			std::string description;
			std::string argName;
			bool optional;
			size_t required;

			Option(char shortName,
				   const std::string& longName = "",
				   const std::string& description = "",
				   const std::string& argName = "",
				   bool optional = false,
				   size_t required = 0)
				: shortName(shortName), longName(longName),
				  description(description), argName(argName),
				  optional(optional), required(required) {
			}

			virtual ~Option() {}

			// Returns true if this option takes one optional argument.
			virtual bool takesOptionalArg() const { return optional; }
		
			// Returns the number of arguments required by this option.
			virtual size_t numRequiredArgs() const { return required; }

			// Invoke this option with arguments given by ARGS
			virtual void process(std::vector<std::string>& args) { return; }

			// Get the description string for this option
			virtual std::string getDescription() const { return description; }
		
			// Print the usage information for this option to STRM
			virtual std::string getHelp(size_t argColWidth = 30,
										size_t lineWidth = 80) const;
		
			std::string getName() const;
		};
 
		template<typename ArgType>
		class StoreOption : public Option {
		private:
			ArgType& arg;
			ArgType val;
		public:
			StoreOption(char shortName,
						const std::string& longName,
						const std::string& description,
						ArgType& arg,
						const std::string& argName = "",
						bool optional = false,
						ArgType val = ArgType())
				: Option(shortName, longName, description, argName, optional,
						 optional ? 0 : 1),
				  arg(arg), val(val) {
			}

			virtual std::string getDescription() const {
				return description + " (default: " + util::string::toString(arg) + ")";
			}
		
			virtual void process(std::vector<std::string>& args) {
				if (!args.empty()) {
					try {
						arg = boost::lexical_cast<ArgType, std::string>(args[0]);
					} catch (boost::bad_lexical_cast& e) {
						throw OptionError("Bad argument for option " + getName());
					}
				} else {
					arg = val;
				}
			}
		};

		template<typename ArgType>
		class StoreConstOption : public Option {
		private:
			ArgType& arg;
			ArgType val;
		public:
			StoreConstOption(char shortName,
							 const std::string& longName,
							 const std::string& description,
							 ArgType& arg,
							 ArgType val)
				: Option(shortName, longName, description),
				  arg(arg), val(val) {
			}
		
			virtual void process(std::vector<std::string>& args) {
				arg = val;
			}
		};

		template<typename Container>
		class AppendConstOption : public Option {
			typedef typename Container::value_type ValueType;
		private:
			Container& c;
			ValueType val;
		public:
			AppendConstOption(char shortName,
							  const std::string& longName,
							  const std::string& description,
							  Container& c,
							  ValueType val)
				: Option(shortName, longName, description),
				  c(c), val(val) {
			}
			
			virtual void process(std::vector<std::string>& args) {
				c.push_back(val);
			}
		};
		
		template<typename Container>
		class AppendOption : public Option {
			Container& c;
		public:
			AppendOption(char shortName,
						 const std::string& longName,
						 const std::string& description,
						 Container& c,
						 const std::string& argName = "")
				: Option(shortName, longName, description, argName, false, 1),
				  c(c) {
			}
		
			virtual void process(std::vector<std::string>& args) {
				typedef typename Container::value_type ValueType;
				try {
					c.push_back(boost::lexical_cast<ValueType, std::string>(args[0]));
				} catch (boost::bad_lexical_cast& e) {
					throw OptionError("Bad argument for option " + getName());
				}
			}
		};

		template<typename Callback>
		class CallbackOption : public Option {
			Callback processor;
		public:
			CallbackOption(char shortName,
						   const std::string& longName,
						   const std::string& description,
						   Callback processor,
						   const std::string& argName,
						   bool optional = false,
						   size_t required = 0)
				: Option(shortName, longName, description, argName, optional, required),
				  processor(processor) {
			}

			virtual void process(std::vector<std::string>& args) {
				processor(args);
			}
		};
					  
		template<typename CounterType>
		class CounterOption : public Option {
			CounterType& counter;
		public:
			CounterOption(char shortName,
						  const std::string& longName,
						  const std::string& description,
						  CounterType& counter)
				: Option(shortName, longName, description),
				  counter(counter) {
			}

			virtual void process(std::vector<std::string>& args) {
				++counter;
			}
		};

		class HelpOption : public Option {
		public:
			HelpOption(char shortName, const std::string& longName)
				: Option(shortName, longName, "Display this help message") {
			}
			void process(std::vector<std::string>& args) {
				throw HelpException();
			}
		};
	
		class Parser {
		private:
			std::string usage;
			std::string progDescription;
			std::string progName;
	
			std::vector<Option*> options;
			std::vector<Argument*> arguments;

			static bool isEndOfOptionsFlag(const std::string& str);
			static bool isShortOption(const std::string& str);
			static bool isLongOption(const std::string& str);

			// Returns the option whose short name is SHORTNAME
			Option* findOption(char shortName) const;

			// Returns the option whose long name begins with LONGNAME.
			// Returns NULL if no option begins with LONGNAME or if no
			// option has long name equal to LONGNAME and multiple options
			// begin with LONGNAME
			Option* findOption(const std::string& longName) const;
		
			void collectOptionArgs(Option* opt,
								   std::vector<std::string>& optArgs,
								   const char**& begin, const char**& end) const;
		
			void parseShortOption(const std::string& arg,
								  const char**& begin, const char**& end) const;
		
			void parseLongOption(const std::string& arg,
								 const char**& begin, const char**& end) const;
		
		public:
			Parser(const std::string& usage,
				   const std::string& progDescription);
		
			void parse(const char** begin, const char** end,
					   bool handleErrors = true,
					   std::ostream& errorStream = std::cerr);

			void addOpt(Option* opt) { options.push_back(opt); }
			void addArg(Argument* arg) { arguments.push_back(arg); }
			void printUsage(std::ostream& strm) const;
			void printHelp(std::ostream& strm) const;
		
			template<typename ArgType>
			void addStoreOpt(char shortName,
							 const std::string& longName,
							 const std::string& description,
							 ArgType& arg,
							 const std::string& argName = "",
							 bool optional = false,
							 ArgType val = ArgType());
		
			template<typename ArgType>
			void addStoreConstOpt(char shortName,
								  const std::string& longName,
								  const std::string& description,
								  ArgType& arg,
								  ArgType val);

			void addStoreTrueOpt(char shortName,
								 const std::string& longName,
								 const std::string& description,
								 bool& arg);

			void addStoreFalseOpt(char shortName,
								  const std::string& longName,
								  const std::string& description,
								  bool& arg);

			template<typename Container>
			void addAppendConstOpt(char shortName,
								   const std::string& longName,
								   const std::string& description,
								   Container& c,
								   typename Container::value_type val);
			
			template<typename Container>
			void addAppendOpt(char shortName,
							  const std::string& longName,
							  const std::string& description,
							  Container& c,		 
							  const std::string& argName = "");

			template<typename CounterType>
			void addCounterOpt(char shortName,
							   const std::string& longName,
							   const std::string& description,
							   CounterType& counter);
		
			template<typename Callback>
			void addCallbackOpt(char shortName,
								const std::string& longName,
								const std::string& description,
								Callback processor,
								const std::string& argName = "",
								bool optional = false,
								size_t required = 0);
		
			template<typename ArgType>
			void addStoreArg(const std::string& argName,
							 const std::string& description,
							 ArgType& arg);

			template<typename Collection>
			void addAppendArg(const std::string& argName,
							  const std::string& description,
							  Collection& c,
							  size_t minArgs = 0,
							  size_t maxArgs = std::numeric_limits<size_t>::max());
		
		};


		template<typename ArgType>
		void Parser::addStoreOpt(char shortName,
								 const std::string& longName,
								 const std::string& description,
								 ArgType& arg,
								 const std::string& argName,
								 bool optional,
								 ArgType val) {
			addOpt(new StoreOption<ArgType>(shortName,
											longName,
											description,
											arg,
											argName,
											optional,
											val));
		}
		

		template<typename ArgType>
		void Parser::addStoreConstOpt(char shortName,
									  const std::string& longName,
									  const std::string& description,
									  ArgType& arg,
									  ArgType val) {
			addOpt(new StoreConstOption<ArgType>(shortName,
												 longName,
												 description,
												 arg,
												 val));
		}

		inline void Parser::addStoreTrueOpt(char shortName,
											const std::string& longName,
											const std::string& description,
											bool& arg) {
			addStoreConstOpt(shortName, longName, description, arg, true);
		}

		inline void Parser::addStoreFalseOpt(char shortName,
											 const std::string& longName,
											 const std::string& description,
											 bool& arg) {
			addStoreConstOpt(shortName, longName, description, arg, false);
		}

		template<typename Container>
		void Parser::addAppendConstOpt(char shortName,
									   const std::string& longName,
									   const std::string& description,
									   Container& c,
									   typename Container::value_type val) {
			addOpt(new AppendConstOption<Container>(shortName,
													longName,
													description,
													c,
													val));
		}

		template<typename Container>
		void Parser::addAppendOpt(char shortName,
								  const std::string& longName,
								  const std::string& description,
								  Container& c,		 
								  const std::string& argName) {
			addOpt(new AppendOption<Container>(shortName,
											   longName,
											   description,
											   c,
											   argName));
		}

		template<typename CounterType>
		void Parser::addCounterOpt(char shortName,
								   const std::string& longName,
								   const std::string& description,
								   CounterType& counter) {
			addOpt(new CounterOption<CounterType>(shortName,
												  longName,
												  description,
												  counter));
		}

		template<typename Callback>
		void Parser::addCallbackOpt(char shortName,
									const std::string& longName,
									const std::string& description,
									Callback processor,
									const std::string& argName,
									bool optional,
									size_t required) {
			addOpt(new CallbackOption<Callback>(shortName,
												longName,
												description,
												processor,
												argName,
												optional,
												required));
		}

		template<typename ArgType>
		void Parser::addStoreArg(const std::string& argName,
								 const std::string& description,
								 ArgType& arg) {
			addArg(new StoreArgument<ArgType>(argName, description, arg));
		}

		template<typename Collection>
		void Parser::addAppendArg(const std::string& argName,
								  const std::string& description,
								  Collection& c,
								  size_t minArgs,
								  size_t maxArgs) {
			addArg(new AppendArgument<Collection>(argName, description, c,
												  minArgs, maxArgs));
		}
	};	
};

#endif // __OPTIONS_HH__
