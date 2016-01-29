#! /usr/bin/env python3
# core.py
# Gao Wang (c) 2015
import os
import numpy as np
import pandas as pd
import deepdish as dd
from .utils_io import InputChecker, map2pandas, load_ddm
from .utils import is_empty, env
from .pyeqtlbma import BFCalculator
from .mix_opt import mixIP
from .utils_workhorse import PosteriorController

def test_association(params):
    params = InputChecker('test_association').apply(params)
    try:
        fn = params["vectors"]["sumstats"]
        sumstats = load_ddm(fn[0], fn[1])
        env.log("Use existing summary statistics data from [{}]".\
                format(os.path.join(os.path.splitext(fn[0])[0], fn[1])))
    except:
        sumstats = {'':{'':{'':{'':0}}}}
    exe = BFCalculator(params["string"], params["int"], params["float"], params["vectors"])
    exe.apply(sumstats, dd.io.load(params["string"]["priors"]))
    res = {"SumStats":
           map2pandas(exe.GetSstats(), "ddm",
                      rownames = exe.GetSstatsRownames(),
                      colnames = ("maf", "n", "pve", "sigmahat", "betahat.geno",
                                  "sebetahat.geno", "betapval.geno")),
           "Abfs":
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
    dd.io.save(params["string"]["output"],
               dict((k, v) for k, v in res.items() if not is_empty(v)),
               compression=("zlib", 9))

def fit_hm(params):
    params = InputChecker('fit_hm').apply(params)
    data = pd.concat(dd.io.load(params["output"], '/' + params["table_abf"])).\
      rename(columns = {'nb_groups' : 'null'})
    data['null'] = 0
    # FIXME: need to implement penalized null
    # FIXME: need to use BF at original scale
    res, converged = mixIP(np.power(10, data), control = params["optimizer_control"])
    if not converged:
        env.error("Convex optimization for hierarchical Model did not converge!")
    # FIXME: eventually all output should be in the same file after I update deepdish with append mode
    dd.io.save(params["output_2"], {params["table_pi"]: pd.Series(res, index = data.columns)},
               compression=("zlib", 9))

def calculate_posterior(params):
    params = InputChecker('calculate_posterior').apply(params)
    # FIXME: all these names
    pc = PosteriorController((params["output"], "/JoinSstats"), (params["output"], "/Abfs"),
                             (params["priors"], "/"), (params["output_2"], "/" + params["table_pi"]),
                             (params["output_3"], "/"), params)
    pc.ScanBlocks()
