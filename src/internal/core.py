#! /usr/bin/env python3
# core.py
# Gao Wang (c) 2015
from .utils_io import *
from .utils import is_empty
from .pyeqtlbma import *

def test_association(params):
    params = check_input(params)
    params = convert_data(params)
    exe = eQtlBma()
    exe.eqtlbma_bf(params["string"], params["int"], params["float"], params["vectors"],
                   {'':{'':{'':{'':0}}}} if "sumstats" in params["None"]
                   else load_sumstats(params["vectors"]["sumstats"][0], params["vectors"]["sumstats"][1]),
                   dd.io.load(params["string"]["priors"]))
    res = {"SumStats":
           convert_data(exe.GetSstats(), "ddm",
                         rownames = exe.GetSstatsRownames(),
                         colnames = ("maf", "n", "pve", "sigmahat", "betahat.geno",
                                     "sebetahat.geno", "betapval.geno")),
           "Abfs":
           convert_data(exe.GetAbfs(), "dm",
                         rownames = exe.GetAbfsNames()),
           "SepPermPvals":
           convert_data(exe.GetSepPermPvals(), "dm",
                         rownames = exe.GetSepPermPvalsRownames()),
           "JoinPermPvals":
           convert_data(exe.GetJoinPermPvals(), "m",
                         rownames = exe.GetJoinPermPvalsRownames())
           }
    dd.io.save(params["string"]["output"],
               dict((k, v) for k, v in res.items() if not is_empty(v)),
               compression=("zlib", 9))
