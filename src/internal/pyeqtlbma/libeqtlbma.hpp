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
typedef std::map<std::string, std::map<std::string, std::vector<std::vector<double> > > > dict_dict_matrixf;

// copied from eqtlbma_bf.cpp
std::vector<quantgen::Snp*> getSnpsForChr(
    const std::map<std::string, std::vector<quantgen::Snp*> > & mChr2VecPtSnps,
    const std::string & chromosome);

// copied from eqtlbma_bf.cpp
void testForAssociations(
    const bool & hasDataNotSstats,
    const std::map<std::string,std::vector<quantgen::Snp*> > & mChr2VecPtSnps,
    const std::string & anchor,
    const size_t & radius,
    const std::vector<std::string> & subgroups,
    const quantgen::Samples & samples,
    const std::string & likelihood,
    const std::string & analysis,
    const bool & need_qnorm,
    const quantgen::Covariates & covariates,
    const quantgen::Grid & iGridL,
    const quantgen::Grid & iGridS,
    const quantgen::PriorMatrices & iPriorM,
    const std::string & bfs,
    const std::string & error_model,
    const float & prop_cov_errors,
    const int & verbose,
    std::map<std::string, quantgen::Gene>::iterator & itG_begin,
    std::map<std::string, quantgen::Gene>::iterator & itG_end,
    size_t & nbAnalyzedGenes,
    size_t & nbAnalyzedPairs);

// copied from eqtlbma_bf.cpp
void makePermutationsSep(
    const std::vector<std::string> & subgroups,
    const quantgen::Samples & samples,
    const std::string & likelihood,
    const bool & need_qnorm,
    const quantgen::Covariates & covariates,
    size_t & nb_permutations,
    const size_t & seed,
    const int & trick,
    const size_t & trick_cutoff,
    const int & perm_sep,
    const gsl_rng * rngPerm,
    const gsl_rng * rngTrick,
    std::map<std::string, quantgen::Gene>::iterator & itG_begin,
    std::map<std::string, quantgen::Gene>::iterator & itG_end);


// copied from eqtlbma_bf.cpp
void makePermutationsJoin(
    const std::vector<std::string> & subgroups,
    const quantgen::Samples & samples,
    const std::string & likelihood,
    const bool & need_qnorm,
    const quantgen::Covariates & covariates,
    const quantgen::Grid & iGridL,
    const quantgen::Grid & iGridS,
    const quantgen::PriorMatrices & iPriorM,
    const std::string & error_model,
    const float & prop_cov_errors,
    size_t & nb_permutations,
    const size_t & seed,
    const int & trick,
    const size_t & trick_cutoff,
    const std::string & permbf,
    const bool & use_max_bf,
    const gsl_rng * rngPerm,
    const gsl_rng * rngTrick,
    std::map<std::string, quantgen::Gene>::iterator & itG_begin,
    std::map<std::string, quantgen::Gene>::iterator & itG_end);

// copied from eqtlbma_bf.cpp
void makePermutations(
    const std::vector<std::string> & subgroups,
    const quantgen::Samples & samples,
    const std::string & likelihood,
    const std::string & analysis,
    const bool & need_qnorm,
    const quantgen::Covariates & covariates,
    const quantgen::Grid & iGridL,
    const quantgen::Grid & iGridS,
    const quantgen::PriorMatrices & iPriorM,
    const std::string & error_model,
    const float & prop_cov_errors,
    size_t & nb_permutations,
    const size_t & seed,
    const int & trick,
    const size_t & trick_cutoff,
    const int & perm_sep,
    const std::string & permbf,
    const bool & use_max_bf,
    std::map<std::string, quantgen::Gene>::iterator & itG_begin,
    std::map<std::string, quantgen::Gene>::iterator & itG_end);

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
