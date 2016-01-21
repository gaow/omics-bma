#! /usr/bin/env python3
# core.py
# Gao Wang (c) 2015
import pandas as pd
import deepdish as dd
from .utils_io import InputChecker, map2pandas, load_ddm
from .utils import is_empty
from .pyeqtlbma import BFCalculator
from .mix_opt import mixIP
from .utils_math import PosteriorCalculator

def test_association(params):
    params = InputChecker('test_association').apply(params)
    sumstats = {'':{'':{'':{'':0}}}}
    if "sumstats" not in params["None"]:
        sumstats = load_ddm(params["vectors"]["sumstats"][0], params["vectors"]["sumstats"][1])
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
    data = dd.io.load(params["output"], params["table_abf"])

def posterior_inference(param):
    pass
