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
#include <gsl/gsl_combination.h>
#include <omp.h>

#include "utils/utils_io.hpp"
#include "utils/utils_math.hpp"
#include "quantgen/data_loader.hpp"

#include "quantgen/gene.hpp"
#include "quantgen/snp.hpp"
#include "quantgen/gene_snp_pair.hpp"

namespace pyeqtlbma {

// copied from eqtlbma_bf.cpp
std::vector<quantgen::Snp *> getSnpsForChr(
	const std::map<std::string, std::vector<quantgen::Snp *> > & mChr2VecPtSnps,
	const std::string & chromosome);

// copied from eqtlbma_bf.cpp
void testForAssociations(
	const bool & hasDataNotSstats,
	const std::map<std::string, std::vector<quantgen::Snp *> > & mChr2VecPtSnps,
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
	const std::vector<std::string> & bfs,
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
	const size_t & nb_permutations,
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
	const size_t & nb_permutations,
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
	const size_t & nb_permutations,
	const size_t & seed,
	const int & trick,
	const size_t & trick_cutoff,
	const int & perm_sep,
	const std::string & permbf,
	const bool & use_max_bf,
	std::map<std::string, quantgen::Gene>::iterator & itG_begin,
	std::map<std::string, quantgen::Gene>::iterator & itG_end);

void extractResSstats(
	const std::vector<std::string> & subgroups,
	const std::map<std::string, quantgen::Gene>::iterator & itG_begin,
	const std::map<std::string, quantgen::Gene>::iterator & itG_end,
	const std::map<std::string, quantgen::Snp> & snp2object,
	std::map<std::string,
	         std::map<std::string,
	                  std::vector<std::vector<double> > > > & res_data,
	std::map<std::string,
	         std::map<std::string,
	                  std::vector<std::string > > > & res_snps);


void extractResSepPermPvalMultiGroup(
	const std::map<std::string, quantgen::Gene>::iterator & itG_begin,
	const std::map<std::string, quantgen::Gene>::iterator & itG_end,
	const std::vector<std::string> & subgroups,
	std::map<std::string, std::vector<std::vector<double> > > & res_data,
	std::map<std::string, std::vector<std::string> > & res_sbgrps
    );

void extractResSepPermPvalSingleGroup(
	const std::map<std::string, quantgen::Gene>::iterator & itG_begin,
	const std::map<std::string, quantgen::Gene>::iterator & itG_end,
	std::map<std::string, std::vector<std::vector<double> > > & res_data
    );

void extractResAbfs(
	const std::map<std::string, quantgen::Gene>::iterator & itG_begin,
	const std::map<std::string, quantgen::Gene>::iterator & itG_end,
	const size_t & nb_subgroups,
	const quantgen::Grid & iGridL,
	const quantgen::Grid & iGridS,
	const quantgen::PriorMatrices & iPriorM,
	const std::vector<std::string> & bfs,
	std::map<std::string, std::vector<std::vector<double> > > & res_data,
	std::map<std::string, std::vector<std::string> > & res_names
    );

void extractResJoinPermPval(
	const std::map<std::string, quantgen::Gene>::iterator & itG_begin,
	const std::map<std::string, quantgen::Gene>::iterator & itG_end,
	const std::string & permbf,
	const bool & use_max_bf,
	std::vector<std::vector<double> > & res_data,
	std::vector<std::string> & res_names
    );

void formatSummaryStats(const std::map<std::string,
	                                   std::map<std::string,
	                                            std::map<std::string,
	                                                     std::map<std::string,
	                                                              double> > > > & sstats,
	const int & verbose,
	std::map<std::string, quantgen::Gene> & gene2object,
	std::map<std::string, quantgen::Snp> & snp2object);

}

#endif
