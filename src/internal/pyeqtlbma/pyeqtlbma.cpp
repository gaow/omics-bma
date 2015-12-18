// pyeqtlbma.cpp
// Gao Wang (c) 2015
#include <iomanip>
#include "pyeqtlbma.hpp"
using namespace std;
using namespace quantgen;
using namespace utils;

namespace pyeqtlbma {
void initialize()
{
}


// implementation follows from the logic in eqtlbma_bf.cpp
int eQtlBma::eqtlbma_bf(
                        const dict_string & param_s,
                        const dict_int & param_i,
                        const dict_float & param_f,
                        const dict_dict_matrixf & sstats,
                        const dict_vectors & param_vs
                        )
{
	set<string> sSnpsToKeep;
	if (param_s.find("snp") != param_s.end()) {
		loadSnpsToKeep(param_s.at("snp"), param_i.at("verbose"), sSnpsToKeep);
		if (sSnpsToKeep.empty())
			return 0;
	}

	vector<string> subgroups;
	Samples samples;
	map<string, Snp> snp2object;
	map<string, vector<Snp *> > mChr2VecPtSnps;
	Covariates covariates;
	map<string, Gene> gene2object;

	// FIXME if(sstats.empty())
	if (1)
		loadRawInputData(param_s.at("geno"), param_s.at("scoord"),
			param_s.at("exp"),
			param_s.at("gcoord"), param_s.at("anchor"), (size_t)param_i.at(
				"cis"),
			param_f.at("maf"), param_s.at("covar"),
			param_s.at("error"), param_vs.at("sbgrp"), sSnpsToKeep,
			param_i.at("verbose"),
			subgroups, samples, snp2object, mChr2VecPtSnps,
			covariates, gene2object);
	else
		// loadSummaryStats(sstats, param_i.at("verbose"), subgroups, gene2object,
		//  snp2object);
		// FIXME: summary statistics have to be provided as map object not files
		return 0;
	if (gene2object.empty() || snp2object.empty())
		return 0;

	// load grid / customized priors
	Grid iGridL(param_s.at("gridL"), true, param_i.at("verbose"));
	Grid iGridS(param_s.at("gridS"), false, param_i.at("verbose"));
	PriorMatrices iPriorM(param_s.at("priorM"), param_s.at("gridM"),
	                      subgroups.size(),
	                      param_i.at("verbose"));
	if (iPriorM.Wg_scalars.size() > iGridL.size()) {
		cerr	<<
		    "WARNING: gridM is larger than gridL; thus will be truncated to the size of gridL"
		        << endl;
	}

	size_t countGenes = 0;
	size_t totalGenes = gene2object.size();
	bool is_perm = param_i.at("nperm") > 0 &&
	               (param_i.at("permsep") != 0 || param_s.at("pbf") != "none");
	if (param_i.at("verbose") > 0) {
		cout	<< "test for association between each pair gene-SNP ..." << endl
		        << "analysis=" << param_s.at("analys")
		        << " likelihood=" << param_s.at("lik")
		        << " error_model=" << param_s.at("error");
		if (param_s.at("error") != "uvlr") // i.e. if 'mvlr' or 'hybrid'
			cout << " prop_cov_errors=" << param_f.at("fiterr");
		// FIXME sstats.empty()
		if (1)
			cout << " anchor=" << param_s.at("anchor") << " radius=" <<
			    param_i.at("cis");
		if (is_perm) {
			cout	<< endl << "permutation" <<
			(param_i.at("nperm") > 1 ? "s=" : "=")
			        << param_i.at("nperm") << " seed=" << param_i.at("seed");
			if (param_i.at("trick") != 0) {
				cout	<< " trick=" << param_i.at("trick")
				        << " trick_cutoff=" << param_i.at("tricut");
			}
			if (param_s.at("analys") == "sep")
				cout << " perm_sep=" << param_i.at("permsep");
			else if (param_s.at("analys") == "join")
				cout << " perm_bf=" << param_s.at("pbf");
			cout << " threads=" << param_i.at("thread");
		}
		cout << endl << flush;
	}
	clock_t startTime = clock();
	size_t nbAnalyzedGenes = 0, nbAnalyzedPairs = 0;
	if (param_i.at("verbose") == 1)
		progressBar("", 0, totalGenes);
	for (map<string, Gene>::iterator itG = gene2object.begin();
	     itG != gene2object.end(); ) {
		map<string, Gene>::iterator itG_begin = itG;
		size_t step_size = min((int)distance(itG, gene2object.end()),
			param_i.at("wrtsize"));
		advance(itG, step_size);
		// FIXME sstats.empty()
		testForAssociations(1, mChr2VecPtSnps, param_s.at(
				"anchor"),
			(size_t)param_i.at("cis"),
			subgroups, samples, param_s.at("lik"), param_s.at("analys"),
			(bool)param_i.at(
				"qnorm"), covariates, iGridL, iGridS, iPriorM,
			param_s.at("bfs"),
			param_s.at("error"), param_f.at("fiterr"), param_i.at("verbose"),
			itG_begin, itG, nbAnalyzedGenes, nbAnalyzedPairs);
		if (is_perm) {
			omp_set_num_threads(param_i.at("thread"));
			makePermutations(subgroups, samples, param_s.at("lik"),
				param_s.at("analys"),
				(bool)param_i.at(
					"qnorm"), covariates, iGridL, iGridS, iPriorM,
				param_s.at("error"), param_f.at("fiterr"),
				(size_t)param_i.at("nperm"), (size_t)param_i.at("seed"),
				param_i.at("trick"), (size_t)param_i.at("tricut"),
				param_i.at("permsep"), param_s.at("pbf"), (bool)param_i.at(
					"maxbf"), itG_begin,
				itG);
		}

		// write results
		// FIXME: sstat empty
		if (param_s.at("analys") == "sep" ||
		    (param_s.at("analys") == "join" && 1 &&
		     param_s.at("error") != "mvlr"))
			extractResSstats(subgroups, itG_begin, itG, snp2object, m_sstats,
				m_sstats_rownames);

		if (param_s.at("analys") == "sep" && param_i.at("nperm") > 0 &&
		    param_i.at("permsep") != 0) {
			if (param_i.at("permsep") == 1)
				extractResSepPermPvalSingleGroup(itG_begin, itG, (size_t)param_i.at(
						"seed"), m_sep_perm_pvals);
			else if (param_i.at("permsep") == 2)
				extractResSepPermPvalMultiGroup(itG_begin, itG, subgroups,
					(size_t)param_i.at(
						"seed"), m_sep_perm_pvals, m_sep_perm_pvals_rownames);
		}


		if (param_s.at("analys") == "join") {
			extractResAbfs(itG_begin, itG, subgroups.size(),
				iGridL, iGridS, iPriorM, param_s.at("bfs"), m_abfs,
				m_abfs_names);
		}

		if (param_s.at("analys") == "join" && param_i.at("nperm") > 0)
			extractResJoinPermPval(itG_begin, itG, (size_t)param_i.at(
					"seed"),
				param_s.at("pbf"),
				(bool)param_i.at(
					"maxbf"), m_join_perm_pvals,
				m_join_perm_pvals_rownames);

		// progress tracker
		countGenes += step_size;
		if (param_i.at("verbose") == 1)
			progressBar("", countGenes, totalGenes);
		// delete processed object to release RAM
		gene2object.erase(itG_begin, itG);
	}
	if (param_i.at("verbose") > 0) {
		if (param_i.at("verbose") == 1)
			cout	<< " (" << fixed << setprecision(2) << getElapsedTime(
				startTime)
			        << " sec)" << endl << flush;
		cout	<< "nb of analyzed gene-SNP pairs: " << nbAnalyzedPairs
		        << " (" << nbAnalyzedGenes << " genes)" << endl;

		cout	<< "END " << __func__
		        << "max.mem -> " << getMaxMemUsedByProcess2Str() << endl;
	}
	return 0;
}


}
