// libqtlbma.hpp
// Gao Wang (c) 2015

#ifndef __LIBEQTLBMA_H_
#define __LIBEQTLBMA_H_

#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <sys/stat.h>

#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <limits>
#include <sstream>
#include <numeric>

#include <gsl/gsl_rng.h>

#include "utils/utils_io.hpp"
#include "utils/utils_math.hpp"

#include "quantgen/gene.hpp"
#include "quantgen/snp.hpp"
#include "quantgen/gene_snp_pair.hpp"

namespace pyeqtlbma {

int eqtlbma_bf();

/// exception handler. Exceptions will be passed to Python.
class Exception
{
public:
	/// constructor
	/// \param msg error message
	Exception(const std::string & msg) : m_msg(msg)
	{
	}


	/// return error message
	const char * message()
	{
		return m_msg.c_str();
	}


	virtual ~Exception()
	{
	};

private:
	/// error message
	std::string m_msg;
};

/// exception, thrown if out of memory
class StopIteration : public Exception
{
public:
	StopIteration(const std::string msg) : Exception(msg)
	{
	};
};


/// exception, thrown if index out of range
class IndexError : public Exception
{
public:
	IndexError(const std::string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if value of range etc
class ValueError : public Exception
{
public:
	ValueError(const std::string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if system error occurs
class SystemError : public Exception
{
public:
	SystemError(const std::string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if a runtime error occurs
class RuntimeError : public Exception
{
public:
	RuntimeError(const std::string msg) : Exception(msg)
	{
	};
};

// initialize C++ module, currently does nothing
void initialize();

}

#endif
