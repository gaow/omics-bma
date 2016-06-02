#! /usr/bin/env python3
__author__ = "Gao Wang"
__copyright__ = "Copyright 2016, Stephens lab"
__email__ = "gaow@uchicago.edu"
__license__ = "MIT"
__version__ = "0.1.0"
import os
import numpy as np
import pandas as pd
import yaml
from .deepdishio import save, load
from .io import ArgumentLoader, dict2map, map2pandas, load_ddm
from .utils import env, is_empty
from .pyeqtlbma import BFCalculator
from .mix_opt import mixIP
from .posterior import PosteriorController

def test_association(prior_data, sumstats_data, output_file, params,
                     genotype_list = None, phenotype_list = None, covariate_list = None,
                     snp_list = None, block_list = None):
    if isinstance(params, str):
        params = yaml.load(open(params))
    params = ArgumentLoader(params)
    env.verbosity = params['verbose']
    params.LoadInputData(genotype_list, phenotype_list, covariate_list, snp_list, block_list)
    if os.path.isfile(sumstats_data):
        sumstats = load_ddm(sumstats_data, '/')
        env.logger.info("Use existing summary statistics data from ``{}``".\
                format(os.path.splitext(sumstats_data)[0]))
    else:
        sumstats = {'':{'':{'':{'':0}}}}
        params.CheckInputData()
    env.logger.debug("\n" + params.Dump())
    params = dict2map(params)
    calculator = BFCalculator(params["string"], params["int"], params["float"], params["vectors"])
    calculator.apply(sumstats, load(prior_data))
    res = {"log10BFs":
           map2pandas(calculator.GetAbfs(), "dm",
                      rownames = calculator.GetAbfsNames()),
           "SepPermPvals":
           map2pandas(calculator.GetSepPermPvals(), "dm",
                      rownames = calculator.GetSepPermPvalsRownames()),
           "JoinPermPvals":
           map2pandas(calculator.GetJoinPermPvals(), "m",
                      rownames = calculator.GetJoinPermPvalsRownames()),
           "JoinSstats":
           map2pandas(calculator.GetJoinSstats(), "ddm",
                      rownames = calculator.GetJoinSstatsRownames())
           }
    save(output_file, dict((k, v) for k, v in res.items() if not is_empty(v)), compression=("zlib", 9))
    if not os.path.isfile(sumstats_data):
        save(sumstats_data,
             map2pandas(calculator.GetSstats(), "ddm",
                        rownames = calculator.GetSstatsRownames(),
                        colnames = ("maf", "n", "pve", "sigmahat", "betahat.geno",
                                    "sebetahat.geno", "betapval.geno")),
             compression = ("zlib", 9))

def fit_hm(association_data, output_file, params):
    if isinstance(params, str):
        params = yaml.load(open(params))
    params = ArgumentLoader(params)
    data = pd.concat(load(association_data, '/log10BFs')).rename(columns = {'nb_groups' : 'null'})
    # FIXME: need to implement penalized null
    data['null'] = 0
    if params['extract_average_bf_per_class']:
        data = data[[x for x in data.columns if not x.endswith('.avg')]]
    # scale data
    data = data.subtract(data.max(axis = 1), axis = 0)
    res, converged = mixIP(np.power(10, data), control = params["optimizer_control"])
    if not converged:
        env.logger.error("Convex optimization for hierarchical Model did not converge!")
    save(output_file, {"pi": pd.Series(res, index = data.columns)}, compression=("zlib", 9))

def calculate_posterior(prior_data, mixture_data, association_data, output_file, params):
    if isinstance(params, str):
        params = yaml.load(open(params))
    params = ArgumentLoader(params)
    pc = PosteriorController((association_data, "/JoinSstats"),
                             (association_data, "/log10BFs"),
                             (prior_data, "/"),
                             (mixture_data, "/pi"),
                             (output_file, "/"), params)
    pc()
