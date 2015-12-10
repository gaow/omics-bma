// libqtlbma.cpp
// Gao Wang (c) 2015
#include "libeqtlbma.hpp"

using namespace std;
using namespace utils;
using namespace quantgen;

vector<Snp *> pyeqtlbma::getSnpsForChr(
                                       const map<string,
                                                 vector<Snp *> > & mChr2VecPtSnps,
                                       const string & chromosome)
{
	vector<Snp *> vPtSnps;
	map<string, vector<Snp *> >::const_iterator itVecPtSnps =
	    mChr2VecPtSnps.find(chromosome);
	if (itVecPtSnps != mChr2VecPtSnps.end())
		vPtSnps = itVecPtSnps->second;
	return(vPtSnps);
}


void pyeqtlbma::testForAssociations(
                                    const bool & hasDataNotSstats,
                                    const map<string,
                                              vector<Snp *> > & mChr2VecPtSnps,
                                    const string & anchor,
                                    const size_t & radius,
                                    const vector<string> & subgroups,
                                    const Samples & samples,
                                    const string & likelihood,
                                    const string & analysis,
                                    const bool & need_qnorm,
                                    const Covariates & covariates,
                                    const Grid & iGridL,
                                    const Grid & iGridS,
                                    const PriorMatrices & iPriorM,
                                    const string & bfs,
                                    const string & error_model,
                                    const float & prop_cov_errors,
                                    const int & verbose,
                                    map<string, Gene>::iterator & itG_begin,
                                    map<string, Gene>::iterator & itG_end,
                                    size_t & nbAnalyzedGenes,
                                    size_t & nbAnalyzedPairs)
{
	for (map<string, Gene>::iterator itG = itG_begin;
	     itG != itG_end; ++itG) {
		if (verbose > 1)
			cout	<< "gene " << itG->first << " (chr "
			        << itG->second.GetChromosome() << ")" << endl;

		if (hasDataNotSstats) {
			if (verbose > 1) {
				vector<Snp *> vPtSnps = getSnpsForChr(mChr2VecPtSnps,
					itG->second.GetChromosome());
				size_t nbSnps = vPtSnps.size();
				cout	<< "chr " << itG->second.GetChromosome() << ": "
				        << nbSnps << " SNP" << (nbSnps > 0 ? "s" : "") << endl;
			}
			itG->second.SetCisSnps(mChr2VecPtSnps, anchor, radius);
			if (verbose > 1)
				cout	<< "gene '" << itG->first << "': " <<
				    itG->second.GetNbCisSnps()
				        << " cis SNP(s)" << endl;
			if (!itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup()) {
				if (verbose > 1)
					cerr	<< "WARNING: skip gene " << itG->second.GetName()
					        << " because it has no SNP in cis" << endl;
				continue;
			}
		}

		if (analysis == "join" && error_model != "uvlr"
		    && !itG->second.HasExplevelsInAllSubgroups(subgroups)) {
			if (verbose > 1)
				cerr	<< "WARNING: skip gene " << itG->second.GetName()
				        << " because option --error " << error_model
				        << " requires expression levels in all subgroups" <<
				    endl;
			continue;
		}
		itG->second.TestForAssociations(hasDataNotSstats, subgroups, samples,
			likelihood, analysis, need_qnorm,
			covariates, iGridL, iGridS, iPriorM, bfs,
			error_model, prop_cov_errors, verbose - 1);
		++nbAnalyzedGenes;
		nbAnalyzedPairs += itG->second.GetNbGeneSnpPairs();
	}

}


void pyeqtlbma::makePermutationsSep(
                                    const vector<string> & subgroups,
                                    const Samples & samples,
                                    const string & likelihood,
                                    const bool & need_qnorm,
                                    const Covariates & covariates,
                                    const size_t & nb_permutations,
                                    const size_t & seed,
                                    const int & trick,
                                    const size_t & trick_cutoff,
                                    const int & perm_sep,
                                    const gsl_rng * rngPerm,
                                    const gsl_rng * rngTrick,
                                    map<string, Gene>::iterator & itG_begin,
                                    map<string, Gene>::iterator & itG_end)
{
	if (perm_sep == 1) {
		gsl_rng_set(rngPerm, seed);
		if (trick != 0)
			gsl_rng_set(rngTrick, seed);
		for (map<string, Gene>::iterator itG = itG_begin;
		     itG != itG_end; ++itG) {
			if (!itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup())
				continue;
			itG->second.MakePermutationsSepAllSubgroups(subgroups, samples,
				likelihood,
				need_qnorm, covariates,
				nb_permutations, trick,
				trick_cutoff,
				rngPerm, rngTrick);
		}
	}else if (perm_sep == 2) {
		for (vector<string>::const_iterator it_sbgrp = subgroups.begin();
		     it_sbgrp != subgroups.end(); ++it_sbgrp) {
			gsl_rng_set(rngPerm, seed);
			if (trick != 0)
				gsl_rng_set(rngTrick, seed);
			for (map<string, Gene>::iterator itG = itG_begin;
			     itG != itG_end; ++itG) {
				if (!itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup())
					continue;
				itG->second.MakePermutationsSepPerSubgroup(*it_sbgrp, samples,
					likelihood,
					need_qnorm, covariates,
					nb_permutations, trick,
					trick_cutoff,
					rngPerm, rngTrick);
			}
		}
	}
}


void pyeqtlbma::makePermutationsJoin(
                                     const vector<string> & subgroups,
                                     const Samples & samples,
                                     const string & likelihood,
                                     const bool & need_qnorm,
                                     const Covariates & covariates,
                                     const Grid & iGridL,
                                     const Grid & iGridS,
                                     const PriorMatrices & iPriorM,
                                     const string & error_model,
                                     const float & prop_cov_errors,
                                     const size_t & nb_permutations,
                                     const size_t & seed,
                                     const int & trick,
                                     const size_t & trick_cutoff,
                                     const string & permbf,
                                     const bool & use_max_bf,
                                     const gsl_rng * rngPerm,
                                     const gsl_rng * rngTrick,
                                     map<string, Gene>::iterator & itG_begin,
                                     map<string, Gene>::iterator & itG_end)
{
	gsl_rng_set(rngPerm, seed);
	if (trick != 0)
		gsl_rng_set(rngTrick, seed);

	for (map<string, Gene>::iterator itG = itG_begin;
	     itG != itG_end; ++itG) {
		if (!itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup())
			continue;
		if (error_model != "uvlr" &&
		    !itG->second.HasExplevelsInAllSubgroups(subgroups))
			continue;
		itG->second.MakePermutationsJoin(subgroups, samples, likelihood,
			need_qnorm,
			covariates, iGridL, iGridS, iPriorM,
			error_model, prop_cov_errors,
			nb_permutations, trick,
			trick_cutoff, permbf, use_max_bf,
			rngPerm, rngTrick);
	}
}


void pyeqtlbma::makePermutations(
                                 const vector<string> & subgroups,
                                 const Samples & samples,
                                 const string & likelihood,
                                 const string & analysis,
                                 const bool & need_qnorm,
                                 const Covariates & covariates,
                                 const Grid & iGridL,
                                 const Grid & iGridS,
                                 const PriorMatrices & iPriorM,
                                 const string & error_model,
                                 const float & prop_cov_errors,
                                 const size_t & nb_permutations,
                                 const size_t & seed,
                                 const int & trick,
                                 const size_t & trick_cutoff,
                                 const int & perm_sep,
                                 const string & permbf,
                                 const bool & use_max_bf,
                                 map<string, Gene>::iterator & itG_begin,
                                 map<string, Gene>::iterator & itG_end)
{

	gsl_rng * rngPerm = NULL, * rngTrick = NULL;

	gsl_rng_env_setup();
	rngPerm = gsl_rng_alloc(gsl_rng_default);
	if (rngPerm == NULL) {
		cerr << "ERROR: can't allocate memory for the RNG" << endl;
		exit(EXIT_FAILURE);
	}
	if (trick != 0) {
		rngTrick = gsl_rng_alloc(gsl_rng_default);
		if (rngTrick == NULL) {
			cerr << "ERROR: can't allocate memory for the RNG" << endl;
			exit(EXIT_FAILURE);
		}
	}

	if (analysis == "sep" && perm_sep != 0)
		makePermutationsSep(subgroups, samples, likelihood, need_qnorm,
			covariates,
			nb_permutations, seed, trick, trick_cutoff,
			perm_sep, rngPerm, rngTrick, itG_begin, itG_end);
	if (analysis == "join" && permbf != "none")
		makePermutationsJoin(subgroups, samples, likelihood, need_qnorm,
			covariates,
			iGridL, iGridS, iPriorM, error_model, prop_cov_errors,
			nb_permutations, seed, trick, trick_cutoff, permbf,
			use_max_bf, rngPerm, rngTrick, itG_begin, itG_end);

	gsl_rng_free(rngPerm);
	if (trick != 0)
		gsl_rng_free(rngTrick);
}


