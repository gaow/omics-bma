#! /usr/bin/env python3
# core.py
# Gao Wang (c) 2015
import os
import numpy as np
import pandas as pd
import deepdish as dd
from .utils_io import ConfigReader, dict2map, map2pandas, load_ddm
from .utils import is_empty, env
from .pyeqtlbma import BFCalculator
from .mix_opt import mixIP
from .utils_workhorse import PosteriorController

def test_association(params):
    params = dict2map(ConfigReader(params))
    if "input_sumstats_data" not in params['None']:
        sumstats = load_ddm(params["vectors"]["input_sumstats_data"], '/')
        env.log("Use existing summary statistics data from [{}]".\
                format(os.path.splitext(params["vectors"]["input_sumstats_data"])[0]))
    else:
        sumstats = {'':{'':{'':{'':0}}}}
    exe = BFCalculator(params["string"], params["int"], params["float"], params["vectors"])
    exe.apply(sumstats, dd.io.load(params["string"]["prior_data"]))
    res = {"log10BFs":
           map2pandas(exe.GetAbfs(), "dm",
                      rownames = exe.GetAbfsNames()),
           "SepPermPvals":
           map2pandas(exe.GetSepPermPvals(), "dm",
                      rownames = exe.GetSepPermPvalsRownames()),
           "JoinPermPvals":
           map2pandas(exe.GetJoinPermPvals(), "m",
                      rownames = exe.GetJoinPermPvalsRownames()),
           "JoinSstats":
           map2pandas(exe.GetJoinSstats(), "ddm",
                      rownames = exe.GetJoinSstatsRownames())
           }
    dd.io.save(params["string"]["association_data"],
               dict((k, v) for k, v in res.items() if not is_empty(v)),
               compression=("zlib", 9))
    res = map2pandas(exe.GetSstats(), "ddm",
                     rownames = exe.GetSstatsRownames(),
                      colnames = ("maf", "n", "pve", "sigmahat", "betahat.geno",
                                  "sebetahat.geno", "betapval.geno")
                                  )
    dd.io.save(params["string"]["output_sumstats_data"], res, compression = ("zlib", 9))


def fit_hm(params):
    params = ConfigReader(params)
    data = pd.concat(dd.io.load(params["association_data"], '/log10BFs')).\
      rename(columns = {'nb_groups' : 'null'})
    # FIXME: need to implement penalized null
    data['null'] = 0
    if params['extract_average_bf_per_class']:
        data = data[[x for x in data.columns if not x.endswith('.avg')]]
    # down-scale data
    data = data - np.max(data)
    res, converged = mixIP(np.power(10, data), control = params["optimizer_control"])
    if not converged:
        env.error("Convex optimization for hierarchical Model did not converge!")
    dd.io.save(params["association_data"], {'pi': pd.Series(res, index = data.columns)},
               compression=("zlib", 9), mode = 'a')

def calculate_posterior(params):
    params = ConfigReader(params)
    pc = PosteriorController((params["association_data"], "/JoinSstats"),
                             (params["association_data"], "/log10BFs"),
                             (params["prior_data"], "/"),
                             (params["association_data"], "/pi"),
                             (params["posterior_data"], "/"), params)
    pc.ScanBlocks()
