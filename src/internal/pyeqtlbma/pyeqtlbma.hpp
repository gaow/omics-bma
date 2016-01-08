// pyeqtlbma.hpp
// Gao Wang (c) 2015

#ifndef __PYEQTLBMA_H_
#define __PYBEQTLBMA_H_

#include "libeqtlbma.hpp"

namespace pyeqtlbma {

typedef std::vector<int> vectori;
typedef std::vector<std::vector<int> > matrixi;
typedef std::vector<std::string> vectors;
typedef std::vector<double> vectorf;
typedef std::vector<std::vector<double> > matrixf;
typedef std::map<std::string, std::vector<double> > dict_vectorf;
typedef std::map<std::string, std::vector<std::vector<double> > > dict_matrixf;
typedef std::map<std::string, std::vector<int> > dict_vectori;
typedef std::map<std::string, std::vector<std::vector<int> > > dict_matrixi;
typedef std::map<std::string, std::vector<std::string> > dict_vectors;
typedef std::map<std::string, std::string> dict_string;
typedef std::map<std::string, int> dict_int;
typedef std::map<std::string, double> dict_float;
typedef std::map<std::string,
                 std::map<std::string,
                          std::vector<std::vector<double> > > >
    dict_dict_matrixf;
typedef std::map<std::string,
                 std::map<std::string,
                          std::vector<std::string > > >
    dict_dict_vectors;

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

class eQtlBma
{
public:
	eQtlBma () {};
	~eQtlBma() {};

	int eqtlbma_bf(
		const dict_string & param_s,
		const dict_int & param_i,
		const dict_float & param_f,
		const dict_vectors & param_vs,
		const dict_dict_matrixf & sstats,
		const dict_matrixf & priors);

	dict_dict_matrixf GetSstats() { return m_sstats; };
	dict_dict_vectors GetSstatsRownames() { return m_sstats_rownames; };
	dict_matrixf GetSepPermPvals() { return m_sep_perm_pvals; }
	dict_vectors GetSepPermPvalsRownames() { return m_sep_perm_pvals_rownames; }
	matrixf GetJoinPermPvals() { return m_join_perm_pvals; }
	vectors GetJoinPermPvalsRownames(){ return m_join_perm_pvals_rownames; }
	dict_matrixf GetAbfs() { return m_abfs; }
	dict_vectors GetAbfsNames() { return m_abfs_names; }

private:
	dict_dict_matrixf m_sstats;
	dict_dict_vectors m_sstats_rownames;
	dict_matrixf m_sep_perm_pvals;
	dict_vectors m_sep_perm_pvals_rownames;
	matrixf m_join_perm_pvals;
	vectors m_join_perm_pvals_rownames;
	dict_matrixf m_abfs;
	dict_vectors m_abfs_names;
};

}
#endif
