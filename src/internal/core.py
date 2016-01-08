#! /usr/bin/env python3
# core.py
# Gao Wang (c) 2015
from OmicsBMA.utils_io import *
from OmicsBMA.pyeqtlbma import *
import deepdish as dd

def test_association(params):
    params = check_input_(params)
    params = convert_data_(params)
    exe = eQtlBma()
    exe.eqtlbma_bf(params["string"], params["int"], params["float"], params["vectors"],
                   {'':{'':[[]]}} if "sumstats" in params["None"]
                   else load_sumstats_(params["string"]["sumstats"]),
                   dd.io.load(params["string"]["priors"]))
    res = {"SumStats": convert_data_(exe.GetSstats(), "ddm", rownames = exe.GetSstatsRownames()),
           "Abfs": convert_data_(exe.GetAbfs(), "dm", rownames = exe.GetAbfsNames()),
           "SepPermPvals": convert_data_(exe.GetSepPermPvals(), "dm",
                                         rownames = exe.GetSepPermPvalsRownames()),
           "JoinPermPvals": convert_data_(exe.GetJoinPermPvals(), "m",
                                          rownames = exe.GetJoinPermPvalsRownames())
           }
    dd.io.save(params["string"]["output"],
               dict((k, v) for k, v in res.items() if not is_empty_(v)),
               compression=("zlib", 9))
