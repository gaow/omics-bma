// libqtlbma.cpp
// Gao Wang (c) 2015
#include "libeqtlbma.hpp"

using namespace std;
using namespace utils;
using namespace quantgen;

inline bool has_value(const vector<string> & sarray, const string & item)
{
	return find(sarray.begin(), sarray.end(), item) != sarray.end();
}


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
                                    const vector<string> & bfs,
                                    const string & error_model,
                                    const float & prop_cov_errors,
                                    const int & verbose,
                                    map<string, Gene>::iterator & itG_begin,
                                    map<string, Gene>::iterator & itG_end,
                                    map<string,
                                        map<string,
                                            vector<vector<double> > > > & beta_n_cov,
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
		// summary statistics for posterior effect size inference
		if (beta_n_cov.find(itG->second.GetName()) == beta_n_cov.end()) {
			map<string, vector<vector<double> > > tmp;
			beta_n_cov[itG->second.GetName()] = tmp;
		}

		itG->second.CalcSstatsHybrid(subgroups, samples, covariates,
			need_qnorm, prop_cov_errors,
			beta_n_cov[itG->second.GetName()]);
		//
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


void pyeqtlbma::extractResSstats(
                                 const vector<string> & subgroups,
                                 const map<string, Gene>::iterator & itG_begin,
                                 const map<string, Gene>::iterator & itG_end,
                                 const map<string, Snp> & snp2object,
                                 map<string,
                                     map<string,
                                         vector<vector<double> > > > & res_data,
                                 map<string,
                                     map<string,
                                         vector<string > > > & res_snps)
{

	for (map<string, Gene>::const_iterator it_gene = itG_begin;
	     it_gene != itG_end; ++it_gene) {
		if (res_data.find(it_gene->second.GetName()) == res_data.end()) {
			map<string, vector<vector<double> > >  tmp1_;
			map<string, vector<string> >  tmp2_;
			res_data[it_gene->second.GetName()] = tmp1_;
			res_snps[it_gene->second.GetName()] = tmp2_;
		}
		for (size_t s = 0; s < subgroups.size(); ++s) {
			if (!it_gene->second.HasAtLeastOneCisSnp(subgroups[s]))
				continue;
			if (res_data[it_gene->second.GetName()].find(subgroups[s]) ==
			    res_data[it_gene->second.GetName()].end()) {
				res_data[it_gene->second.GetName()][subgroups[s]] = {};
				res_snps[it_gene->second.GetName()][subgroups[s]] = {};
			}
			for (vector<GeneSnpPair>::const_iterator it_pair
			         = it_gene->second.BeginPair();
			     it_pair != it_gene->second.EndPair(); ++it_pair) {
				if (!it_pair->HasResults(subgroups[s]))
					continue;
				res_snps[it_gene->second.GetName()][subgroups[s]].push_back(
					it_pair->GetSnpName());
				vector<double> tmp_ { snp2object.find(it_pair->GetSnpName())->
					                  second.GetMinorAlleleFreq(subgroups[s]),
					                  (double)it_pair->GetSampleSize(
										  subgroups[s]),
					                  it_pair->GetPve(subgroups[s]),
					                  it_pair->GetSigmahat(subgroups[s]),
					                  it_pair->GetBetahatGeno(subgroups[s]),
					                  it_pair->GetSebetahatGeno(subgroups[s]),
					                  it_pair->GetBetapvalGeno(subgroups[s]) };
				res_data[it_gene->second.GetName()][subgroups[s]].push_back(tmp_);
			}   // end of loop over cis snps
		}       // end of loop over subgroups
	}           // end of loop over genes
}


void pyeqtlbma::extractResSepPermPvalMultiGroup(
                                                const map<string,
                                                          Gene>::iterator & itG_begin,
                                                const map<string,
                                                          Gene>::iterator & itG_end,
                                                const vector<string> & subgroups,
                                                map<string,
                                                    vector<vector<double> > > & res_data,
                                                map<string, vector<string> > & res_sbgrps
                                                )
{
	for (map<string, Gene>::const_iterator it_gene = itG_begin;
	     it_gene != itG_end; ++it_gene) {
		for (vector<string>::const_iterator it_sbgrp = subgroups.begin();
		     it_sbgrp != subgroups.end(); ++it_sbgrp) {
			if (it_gene->second.GetNbGeneSnpPairs() > 0) {
				if (res_data.find(it_gene->second.GetName()) ==
				    res_data.end()) {
					res_data[it_gene->second.GetName()] = {};
					res_sbgrps[it_gene->second.GetName()] = {};
				}
				vector<double> tmp_ { (double)it_gene->second.GetNbGeneSnpPairs(
										  *it_sbgrp),
					                  it_gene->second.GetPermutationPvalueSep(
										  *it_sbgrp),
					                  (double)it_gene->second.
					                  GetNbPermutationsSep(
										  *it_sbgrp),
					                  it_gene->second.GetTrueMinPval(*it_sbgrp) };
				res_data[it_gene->second.GetName()].push_back(tmp_);
				res_sbgrps[it_gene->second.GetName()].push_back(*it_sbgrp);
			}
		}
	}
}


void pyeqtlbma::extractResSepPermPvalSingleGroup(
                                                 const map<string,
                                                           Gene>::iterator & itG_begin,
                                                 const map<string,
                                                           Gene>::iterator & itG_end,
                                                 map<string, vector<vector<double> > > & res_data
                                                 )
{
	for (map<string, Gene>::const_iterator it_gene = itG_begin;
	     it_gene != itG_end; ++it_gene) {
		if (it_gene->second.GetNbGeneSnpPairs() > 0) {
			if (res_data.find(it_gene->second.GetName()) == res_data.end()) {
				res_data[it_gene->second.GetName()] = {};
			}
			vector<double> tmp_ { (double)it_gene->second.GetNbGeneSnpPairs(),
				                  it_gene->second.GetPermutationPvalueSep(),
				                  (double)it_gene->second.GetNbPermutationsSep(),
				                  it_gene->second.GetTrueMinPval() };
			res_data[it_gene->second.GetName()].push_back(tmp_);
		}
	}
}


// this function combines writeResAbfsRaw and writeResAbfsAvgGrids
void pyeqtlbma::extractResAbfs(
                               const map<string, Gene>::iterator & itG_begin,
                               const map<string, Gene>::iterator & itG_end,
                               const size_t & nb_subgroups,
                               const Grid & iGridL,
                               const Grid & iGridS,
                               const PriorMatrices & iPriorM,
                               const vector<string> & bfs,
                               const bool & outw,
                               map<string, vector<vector<double> > > & res_data,
                               map<string, vector<string> > & res_names
                               )
{
	stringstream ssConfig;
	gsl_combination * comb;
	bool colnames_saved = false;

	res_names["colnames"] = { "nb_groups" };

	// write results
	for (map<string, Gene>::const_iterator it_gene = itG_begin;
	     it_gene != itG_end; ++it_gene) {
		for (vector<GeneSnpPair>::const_iterator it_pair
		         = it_gene->second.BeginPair();
		     it_pair != it_gene->second.EndPair();
		     ++it_pair) {
			if (res_data.find(it_gene->second.GetName()) ==
			    res_data.end()) {
				res_data[it_gene->second.GetName()] = {};
				res_names[it_gene->second.GetName()] = {};
			}
			res_names[it_gene->second.GetName()].push_back(it_pair->GetSnpName());
			vector<double> tmp_ { (double)it_pair->GetNbSubgroups() };

			// write gen BFs (large grid)
			if (!colnames_saved) {
				for (size_t i = 0; i < iGridL.size(); ++i)
					res_names["colnames"].push_back("gen." + to_string(i + 1));
			}
			for (vector<double>::const_iterator it
			         = it_pair->BeginUnweightedAbf("gen");
			     it != it_pair->EndUnweightedAbf("gen"); ++it)
				tmp_.push_back(*it);

			// write gen-fix BFs (large grid)
			if (!colnames_saved) {
				for (size_t i = 0; i < iGridL.size(); ++i)
					res_names["colnames"].push_back("gen-fix." + to_string(
							i +
							1));
			}
			for (vector<double>::const_iterator it
			         = it_pair->BeginUnweightedAbf("gen-fix");
			     it != it_pair->EndUnweightedAbf("gen-fix"); ++it)
				tmp_.push_back(*it);

			// write gen-maxh BFs (large grid)
			if (!colnames_saved) {
				for (size_t i = 0; i < iGridL.size(); ++i)
					res_names["colnames"].push_back("gen-maxh." +
						to_string(i + 1));
			}
			for (vector<double>::const_iterator it
			         = it_pair->BeginUnweightedAbf("gen-maxh");
			     it != it_pair->EndUnweightedAbf("gen-maxh"); ++it)
				tmp_.push_back(*it);

			// write averaged gen, gen-fix, gen-maxh BFs
			if (outw) {
				if (!colnames_saved) {
					res_names["colnames"].push_back("gen.avg");
					res_names["colnames"].push_back("gen-fix.avg");
					res_names["colnames"].push_back("gen-maxh.avg");
					if (has_value(bfs, "sin") || has_value(bfs, "all"))
						res_names["colnames"].push_back("gen-sin.avg");

					if (has_value(bfs, "all"))
						res_names["colnames"].push_back("all.avg");
					if (has_value(bfs, "customized"))
						res_names["colnames"].push_back("customized.avg");

				}
				tmp_.push_back(it_pair->GetWeightedAbf("gen"));
				tmp_.push_back(it_pair->GetWeightedAbf("gen-fix"));
				tmp_.push_back(it_pair->GetWeightedAbf("gen-maxh"));
				if (has_value(bfs, "sin") || has_value(bfs, "all"))
					tmp_.push_back(it_pair->GetWeightedAbf("gen-sin"));


				if (has_value(bfs,
						"all")) tmp_.push_back(it_pair->GetWeightedAbf("all"));
				if (has_value(bfs, "customized"))
					tmp_.push_back(it_pair->GetWeightedAbf("customized"));
			}

			// write the BFs for each config (small grid)
			if (has_value(bfs, "sin") || has_value(bfs, "all")) {
				for (size_t k = 1; k <= nb_subgroups; ++k) {
					comb = gsl_combination_calloc(nb_subgroups, k);
					if (comb == NULL) {
						cerr	<<
						    "ERROR: can't allocate memory for the combination"
						        << endl;
						exit(EXIT_FAILURE);
					}
					while (true) {
						ssConfig.str("");
						ssConfig << gsl_combination_get(comb, 0) + 1;
						if (comb->k > 1)
							for (size_t i = 1; i < k; ++i)
								ssConfig << "-" <<
								    gsl_combination_get(comb, i) + 1;
						// write individual BFs
						for (size_t j = 0; j < iGridS.size(); ++j) {
							if (!colnames_saved) {
								res_names["colnames"].push_back(
									"cfg." + ssConfig.str() + "." +
									to_string(j + 1));
							}
							tmp_.push_back(*(it_pair->BeginUnweightedAbf(
												 ssConfig.str()) + j));
						}
						// write avg BFs
						if (outw) {
							if (!colnames_saved) res_names["colnames"].push_back(
									"cfg." + ssConfig.str() + ".avg");
							tmp_.push_back(it_pair->GetWeightedAbf(
									ssConfig.str()));
						}
						if (gsl_combination_next(comb) != GSL_SUCCESS)
							break;
					}
					gsl_combination_free(comb);
					if (!has_value(bfs, "all"))
						break;
				}
			}
			// write the BFs for each customized prior (with customized grid)
			if (has_value(bfs, "customized")) {
				for (size_t m = 0; m < iPriorM.Wg_names.size(); ++m) {
					// write indiviual BFs
					for (size_t j = 0; j < iPriorM.Wg_scalars.size(); ++j) {
						if (!colnames_saved) {
							res_names["colnames"].push_back(iPriorM.Wg_names[m] + "." + to_string(
									j + 1));
						}
						tmp_.push_back(*(it_pair->BeginUnweightedAbf(iPriorM.
											 Wg_names[m]) + j));
					}
					// write average BFs
					if (outw) {
						if (!colnames_saved) res_names["colnames"].push_back(
								iPriorM.Wg_names[
								    m] + ".avg");
						tmp_.push_back(it_pair->GetWeightedAbf(iPriorM.Wg_names[
								    m]));
					}
				}
			}
			// collect results data
			res_data[it_gene->second.GetName()].push_back(tmp_);
			if (!colnames_saved) colnames_saved = true;
		}
	}
}


void pyeqtlbma::extractResJoinPermPval(
                                       const map<string,
                                                 Gene>::iterator & itG_begin,
                                       const map<string,
                                                 Gene>::iterator & itG_end,
                                       const string & permbf,
                                       const bool & use_max_bf,
                                       vector<vector<double> > & res_data,
                                       vector<string> & res_names
                                       )
{

	for (map<string, Gene>::const_iterator it_gene = itG_begin;
	     it_gene != itG_end; ++it_gene) {
		if (it_gene->second.GetNbGeneSnpPairs() > 0) {
			res_names.push_back(it_gene->second.GetName());
			vector<double> tmp_ = { (double)it_gene->second.GetNbGeneSnpPairs(),
				                    (double)it_gene->second.
				                    GetPermutationPvalueJoin(),
				                    (double)it_gene->second.
				                    GetNbPermutationsJoin(),
				                    it_gene->second.GetTrueL10Abf(use_max_bf),
				                    it_gene->second.GetMedianPermL10Abf() };
			res_data.push_back(tmp_);
		}
	}
}


void pyeqtlbma::formatSummaryStats(const map<string,
                                             map<string,
                                                 map<string,
                                                     map<string,
                                                         double> > > > & sstats,
                                   const int & verbose,
                                   map<string, Gene> & gene2object,
                                   map<string, Snp> & snp2object)
{

	Gene * pt_gene = NULL;
	Snp * pt_snp = NULL;

	vector<GeneSnpPair>::iterator it_gsp;
	size_t idx_snp = string::npos;

	// loop over genes
	for (auto const & i : sstats) {
		// create gene
		if (gene2object.find(i.first) == gene2object.end())
			gene2object.insert(make_pair(i.first, Gene(i.first)));
		pt_gene = &(gene2object[i.first]);
		// loop over subgroups
		for (auto const & j  : i.second) {
			// loop over SNPs
			for (auto const & k : j.second) {
				// create SNP if necessary
				if (snp2object.find(k.first) == snp2object.end())
					snp2object.insert(make_pair(k.first, Snp(k.first)));
				pt_snp = &(snp2object[k.first]);

				// get gene-snp pair (create it if necessary)
				idx_snp = pt_gene->FindIdxSnp(pt_snp);
				if (idx_snp == string::npos) {
					pt_gene->AddCisSnp(pt_snp);
					it_gsp = pt_gene->AddGeneSnpPair(pt_snp->GetName(), "uvlr");
				} else
					it_gsp = pt_gene->FindGeneSnpPair(idx_snp);

				it_gsp->SetSstats(j.first,
					(size_t)k.second.at("n"),
					k.second.at("sigmahat"),
					k.second.at("betahat.geno"),
					k.second.at("sebetahat.geno"));
			}   // end loop over SNPs
		}       // end loop over subgroups
	}           // end loop over genes
}


