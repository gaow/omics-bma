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
int BFCalculator::apply(const dict_x4_float & sstats,
                        const dict_matrixf & priors)
{
	set<string> sSnpsToKeep;
	if (m_s.find("snp") != m_s.end()) {
		loadSnpsToKeep(m_s.at("snp"), m_i.at("verbose"), sSnpsToKeep);
		if (sSnpsToKeep.empty())
			return 0;
	}

	Samples samples;
	map<string, Snp> snp2object;
	map<string, vector<Snp *> > mChr2VecPtSnps;
	Covariates covariates;
	map<string, Gene> gene2object;

	if (sstats.begin()->first.empty())
		loadRawInputData(m_s.at("geno"), m_s.at("scoord"),
			m_s.at("exp"),
			m_s.at("gcoord"), m_s.at("anchor"), (size_t)m_i.at(
				"cis"),
			m_f.at("maf"), m_s.at("covar"),
			m_s.at("error"), m_vs.at("sbgrp"), sSnpsToKeep,
			m_i.at("verbose"),
			m_subgroups, samples, snp2object, mChr2VecPtSnps,
			covariates, gene2object);
	else
		formatSummaryStats(sstats, m_i.at(
				"verbose"), gene2object,
			snp2object);

	if (gene2object.empty() || snp2object.empty())
		return 0;

	// load grid / customized priors
	Grid iGridL(priors, "gridL", true, m_i.at("verbose"));
	Grid iGridS(priors, "gridS", false, m_i.at("verbose"));
	PriorMatrices iPriorM(priors, "gridM", { "gridL", "gridS" },
	                      m_i.at("verbose"));
	//
	size_t countGenes = 0;
	size_t totalGenes = gene2object.size();
	bool is_perm = m_i.at("nperm") > 0 &&
	               (m_i.at("permsep") != 0 || m_s.at("pbf") != "none");
	if (m_i.at("verbose") > 0) {
		cout	<< "test for association between each pair gene-SNP ..." << endl
		        << "analysis=" << m_s.at("analys")
		        << " likelihood=" << m_s.at("lik")
		        << " error_model=" << m_s.at("error");
		if (m_s.at("error") != "uvlr") // i.e. if 'mvlr' or 'hybrid'
			cout << " prop_cov_errors=" << m_f.at("fiterr");
		if (sstats.begin()->first.empty())
			cout << " anchor=" << m_s.at("anchor") << " radius=" <<
			    m_i.at("cis");
		if (is_perm) {
			cout	<< endl << "permutation" <<
			(m_i.at("nperm") > 1 ? "s=" : "=")
			        << m_i.at("nperm") << " seed=" << m_i.at("seed");
			if (m_i.at("trick") != 0) {
				cout	<< " trick=" << m_i.at("trick")
				        << " trick_cutoff=" << m_i.at("tricut");
			}
			if (m_s.at("analys") == "sep")
				cout << " perm_sep=" << m_i.at("permsep");
			else if (m_s.at("analys") == "join")
				cout << " perm_bf=" << m_s.at("pbf");
			cout << " threads=" << m_i.at("thread");
		}
		cout << endl << flush;
	}
	clock_t startTime = clock();
	size_t nbAnalyzedGenes = 0, nbAnalyzedPairs = 0;
	if (m_i.at("verbose") == 1)
		progressBar("", 0, totalGenes);
	for (map<string, Gene>::iterator itG = gene2object.begin();
	     itG != gene2object.end(); ) {
		map<string, Gene>::iterator itG_begin = itG;
		size_t step_size = min((int)distance(itG, gene2object.end()),
			m_i.at("wrtsize"));
		advance(itG, step_size);
		testForAssociations(
			sstats.begin()->first.empty(), mChr2VecPtSnps, m_s.at(
				"anchor"),
			(size_t)m_i.at("cis"),
			m_subgroups, samples, m_s.at("lik"), m_s.at("analys"),
			(bool)m_i.at(
				"qnorm"), covariates, iGridL, iGridS, iPriorM,
			m_vs.at("bfs"),
			m_s.at("error"), m_f.at("fiterr"), m_i.at("verbose"),
			itG_begin, itG, m_beta_n_cov, nbAnalyzedGenes, nbAnalyzedPairs);
		if (is_perm) {
			omp_set_num_threads(m_i.at("thread"));
			makePermutations(m_subgroups, samples, m_s.at("lik"),
				m_s.at("analys"),
				(bool)m_i.at(
					"qnorm"), covariates, iGridL, iGridS, iPriorM,
				m_s.at("error"), m_f.at("fiterr"),
				(size_t)m_i.at("nperm"), (size_t)m_i.at("seed"),
				m_i.at("trick"), (size_t)m_i.at("tricut"),
				m_i.at("permsep"), m_s.at("pbf"), (bool)m_i.at(
					"maxbf"), itG_begin,
				itG);
		}

		// write results
		if (m_s.at("analys") == "sep" ||
		    (m_s.at("analys") == "join" && sstats.begin()->first.empty() &&
		     m_s.at("error") != "mvlr"))
			extractResSstats(m_subgroups, itG_begin, itG, snp2object, m_sstats,
				m_sstats_rownames);

		if (m_s.at("analys") == "sep" && m_i.at("nperm") > 0 &&
		    m_i.at("permsep") != 0) {
			if (m_i.at("permsep") == 1)
				extractResSepPermPvalSingleGroup(itG_begin, itG,
					m_sep_perm_pvals);
			else if (m_i.at("permsep") == 2)
				extractResSepPermPvalMultiGroup(itG_begin, itG, m_subgroups,
					m_sep_perm_pvals, m_sep_perm_pvals_rownames);
		}


		if (m_s.at("analys") == "join") {
			extractResAbfs(itG_begin, itG, m_subgroups.size(),
				iGridL, iGridS, iPriorM, m_vs.at("bfs"),
				(bool)m_i.at("out_avg"), m_abfs,
				m_abfs_names);
		}

		if (m_s.at("analys") == "join" && m_i.at("nperm") > 0)
			extractResJoinPermPval(itG_begin, itG, m_s.at("pbf"),
				(bool)m_i.at("maxbf"), m_join_perm_pvals,
				m_join_perm_pvals_rownames);

		// progress tracker
		countGenes += step_size;
		if (m_i.at("verbose") == 1)
			progressBar("", countGenes, totalGenes);
		// delete processed object to release RAM
		gene2object.erase(itG_begin, itG);
	}
	if (m_i.at("verbose") > 0) {
		if (m_i.at("verbose") == 1)
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


dict_vectori get_eqtlbma_configurations(const size_t & nb_subgroups,
                                        const bool & is_all)
{
	gsl_combination * comb;
	dict_vectori res;

	for (size_t k = 1; k <= nb_subgroups; ++k) {
		comb = gsl_combination_calloc(nb_subgroups, k);
		if (comb == NULL) {
			cerr	<<
			    "ERROR: can't allocate memory for the combination"
			        << endl;
		}
		while (true) {
			stringstream config_name;
			vectori gamma(nb_subgroups, 0);

			config_name << gsl_combination_get(comb, 0) + 1;
			gamma[gsl_combination_get(comb, 0)] = 1;

			if (comb->k > 1) {
				for (size_t i = 1; i < comb->k; ++i) {
					config_name << "-" << gsl_combination_get(comb, i) + 1;
					gamma[gsl_combination_get(comb, i)] = 1;
				}
			}

			res[config_name.str()] = gamma;

			if (gsl_combination_next(comb) != GSL_SUCCESS)
				break;
		}
		gsl_combination_free(comb);
		if (!is_all)
			break;
	}
	return res;
}


}
