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
from .utils_math import PosteriorCalculator

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
                      rownames = exe.GetJoinPermPvalsRownames())
           }
    dd.io.save(params["string"]["output"],
               dict((k, v) for k, v in res.items() if not is_empty(v)),
               compression=("zlib", 9))

def fit_hm(params):
    params = InputChecker('fit_hm').apply(params)
    data = dd.io.load(params["output"], '/' + params["table_abf"])
    res, converged = mixIP(pd.concat(data), control = params["optimizer_control"])
    if not converged:
        env.error("Convex optimization for hierarchical Model did not converge!")
    # FIXME: eventually all output should be in the same file after I update deepdish with append mode
    dd.io.save(params["output_2"], {params["table_pi"]: np.array(res)}, compression=("zlib", 9))

def posterior_inference(param):
    pass
